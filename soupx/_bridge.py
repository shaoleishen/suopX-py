"""Bridge layer between Python/AnnData and the Rust core engine.

See REFACTOR_PLAN.md §3 for architecture.
"""

from __future__ import annotations

import gzip
import warnings
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from anndata import AnnData

# ---------------------------------------------------------------------------
# Rust core import
# ---------------------------------------------------------------------------
try:
    import soupx_core as _core
    _HAS_RUST_CORE = True
except ImportError as e:
    _HAS_RUST_CORE = False
    _core = None
    warnings.warn(
        f"Rust core (soupx_core) not found: {e}. "
        "Run `maturin develop --release` to build it. "
        "Python fallback not yet available."
    )


# ============================================================================
# Internal helpers
# ============================================================================

def _csr_to_rust_args(mat: sparse.csr_matrix) -> tuple:
    """Convert scipy CSR matrix to (data, indices, indptr, shape) for Rust."""
    if mat.dtype != np.float64:
        mat = mat.astype(np.float64)
    indices = mat.indices.astype(np.int64)
    indptr = mat.indptr.astype(np.int64)
    return (mat.data, indices, indptr, mat.shape)


def _clusters_to_rust(adata: AnnData) -> np.ndarray:
    """Extract cluster labels as integer array."""
    key = adata.uns.get("soupx_cluster_key", "leiden")
    labels = adata.obs[key].values
    # Convert to categorical codes (0, 1, 2, ...)
    if labels.dtype == object or labels.dtype.kind in ('U', 'S'):
        uniq = {v: i for i, v in enumerate(np.unique(labels))}
        return np.array([uniq[v] for v in labels], dtype=np.uint32)
    return labels.astype(np.uint32)


def _ensure_core():
    """Raise if Rust core is not available."""
    if not _HAS_RUST_CORE:
        raise ImportError(
            "Rust core (soupx_core) is required. Run `maturin develop --release`."
        )


# ============================================================================
# load_10x — parse CellRanger output
# ============================================================================

def load_10x(
    adata: AnnData,
    cellranger_dir: str | Path,
    soup_range: tuple[int, int] = (0, 100),
    keep_droplets: bool = False,
) -> None:
    """Load all droplet data from a 10X CellRanger output directory.

    Reads the raw_gene_bc_matrices to get the full droplet matrix,
    estimates the soup profile from empty droplets (UMI in soup_range),
    and stores results in adata.uns and optionally adata.raw.

    Parameters
    ----------
    adata
        AnnData with cell-level counts in adata.X (e.g. from sc.read_10x_mtx).
    cellranger_dir
        Path to CellRanger outs/ directory.
    soup_range
        UMI range for empty droplets (default: 0-100).
    keep_droplets
        If True, retains the full droplet matrix in adata.raw.
        If False (default), stores only the soup profile to save memory.
    """
    cellranger_dir = Path(cellranger_dir)
    raw_dir = _find_raw_matrix_dir(cellranger_dir)

    # Read genes from raw matrix (may differ from filtered)
    features_path = _find_file(raw_dir, "features.tsv", "genes.tsv")
    raw_ids, raw_names = _read_features(features_path)

    # Map raw gene names to AnnData var_names
    adata_genes = list(adata.var_names)
    gene_map = {name: i for i, name in enumerate(raw_names)}
    matched_indices = [gene_map.get(g, -1) for g in adata_genes]
    n_matched = sum(1 for i in matched_indices if i >= 0)
    if n_matched < len(adata_genes):
        warnings.warn(
            f"Only {n_matched}/{len(adata_genes)} genes matched between raw and filtered matrices."
        )

    # Read raw mtx — but only the columns we need
    from scipy.io import mmread
    mtx_path = _find_file(raw_dir, "matrix.mtx")
    if mtx_path.suffix == ".gz":
        with gzip.open(mtx_path, "rt") as f:
            tod_full = sparse.csr_matrix(mmread(f).tocoo())
    else:
        tod_full = sparse.csr_matrix(mmread(str(mtx_path)).tocoo())

    # Transpose if needed (10X format: barcodes × genes → genes × barcodes)
    n_raw_genes, n_raw_barcodes = len(raw_ids), _count_lines(_find_file(raw_dir, "barcodes.tsv"))
    if tod_full.shape[0] == n_raw_barcodes and tod_full.shape[1] == n_raw_genes:
        tod_full = tod_full.T.tocsr()
    tod_full = tod_full.astype(np.float64)

    # Identify empty droplets
    n_umis = np.array(tod_full.sum(axis=0)).ravel()
    empty_mask = (n_umis >= soup_range[0]) & (n_umis <= soup_range[1])
    n_empty = int(np.sum(empty_mask))

    if n_empty == 0:
        warnings.warn(
            f"No empty droplets found in UMI range {soup_range}. Try a wider range."
        )
        empty_counts = np.zeros(len(raw_ids))
    else:
        empty_counts = np.array(tod_full[:, empty_mask].sum(axis=1)).ravel()
    total_empty = empty_counts.sum()

    # Build soup profile aligned to AnnData genes
    full_soup_profile = empty_counts / max(total_empty, 1)
    soup_profile = np.zeros(len(adata_genes))
    for i, raw_idx in enumerate(matched_indices):
        if raw_idx >= 0:
            soup_profile[i] = full_soup_profile[raw_idx]

    _store_soup_profile(adata, soup_profile, empty_counts, soup_range, n_empty)

    if keep_droplets:
        # Subset to matched genes, keep all droplets
        matched_raw = [i for i in matched_indices if i >= 0]
        tod_subset = tod_full[matched_raw, :]
        adata.raw = AnnData(tod_subset.T, dtype=np.float64)
        adata.raw.var_names = adata.var_names
    else:
        adata.uns["soupx_n_droplets"] = tod_full.shape[1]
        adata.uns["soupx_n_empty"] = n_empty

    del tod_full
    gc = __import__("gc")
    gc.collect()


# ============================================================================
# set_clusters
# ============================================================================

def set_clusters(adata: AnnData, key: str = "leiden") -> None:
    """Set the clustering labels for background estimation."""
    if key not in adata.obs:
        raise KeyError(f"Cluster key '{key}' not found in adata.obs")
    adata.uns["soupx_cluster_key"] = key


# ============================================================================
# auto_est_cont — calls Rust core
# ============================================================================

def auto_est_cont(
    adata: AnnData,
    *,
    tfidf_min: float = 0.2,
    soup_quantile: float = 0.90,
    max_markers: int = 100,
    contamination_range: tuple[float, float] = (0.01, 0.8),
    prior_rho: float = 0.05,
    prior_rho_stddev: float = 0.10,
    force_accept: bool = False,
    verbose: bool = True,
) -> None:
    """Automatically estimate contamination fraction ρ via Rust core.

    Stores: adata.obs["soupx_rho"], adata.uns["soupx_fit"].
    """
    _ensure_core()

    if verbose:
        print("auto_est_cont: estimating contamination fraction...")

    # Prepare inputs
    toc = adata.X
    if sparse.issparse(toc):
        toc = toc.tocsr()
    else:
        toc = sparse.csr_matrix(toc)
    data, indices, indptr, shape = _csr_to_rust_args(toc)

    clusters_vec = _clusters_to_rust(adata).tolist()
    soup_profile = adata.uns["soupx_soup_profile"]

    rho = _core.auto_est_cont(
        data, indices, indptr, shape,
        clusters_vec,
        soup_profile.tolist(),
        tfidf_min,
        soup_quantile,
        max_markers,
        contamination_range,
        prior_rho,
        prior_rho_stddev,
    )

    # Store results
    adata.obs["soupx_rho"] = rho
    adata.uns["soupx_fit"] = {
        "rho": rho,
        "prior_rho": prior_rho,
        "prior_rho_stddev": prior_rho_stddev,
        "tfidf_min": tfidf_min,
    }

    if verbose:
        print(f"  Estimated ρ = {rho:.4f}")
        if force_accept:
            print("  (force_accept=True, skipping quality warnings)")


# ============================================================================
# estimate_contamination (manual) — calls Rust core
# ============================================================================

def estimate_contamination(
    adata: AnnData,
    gene_list: dict[str, list[str]],
    *,
    maximum_contamination: float = 1.0,
    fdr: float = 0.05,
    force_accept: bool = False,
) -> float:
    """Manual ρ estimation using user-specified gene sets.

    Parameters
    ----------
    adata
        AnnData with cell counts and soup profile.
    gene_list
        Dict mapping gene set names to lists of gene symbols.
        e.g. {"HB": ["HBB", "HBA2"], "IG": ["IGHG1", "IGKC"]}
    maximum_contamination
        Upper bound for ρ.
    fdr
        False discovery rate for non-expressing cell test.

    Returns
    -------
    Estimated ρ.
    """
    _ensure_core()

    toc = sparse.csr_matrix(adata.X) if not sparse.issparse(adata.X) else adata.X.tocsr()
    data, indices, indptr, shape = _csr_to_rust_args(toc)

    # Map gene names to global indices
    var_names = list(adata.var_names)
    gene_sets: list[list[int]] = []
    for name, genes in gene_list.items():
        indices_list = [var_names.index(g) for g in genes if g in var_names]
        if indices_list:
            gene_sets.append(indices_list)
        else:
            warnings.warn(f"Gene set '{name}': no genes found in adata.var_names")

    if not gene_sets:
        raise ValueError("No valid gene sets found. Check gene names against adata.var_names.")

    clusters_vec = _clusters_to_rust(adata).tolist()
    soup_profile = adata.uns["soupx_soup_profile"]

    rho = _core.estimate_contamination(
        data, indices, indptr, shape,
        gene_sets,
        clusters_vec,
        soup_profile.tolist(),
        maximum_contamination,
        fdr,
    )

    if not force_accept and rho > 0.5:
        warnings.warn(f"Estimated ρ = {rho:.4f} is high. Check gene sets or set force_accept=True.")

    return rho


# ============================================================================
# adjust_counts — calls Rust core
# ============================================================================

def adjust_counts(
    adata: AnnData,
    *,
    method: Literal["subtraction", "soup_only", "multinomial"] = "subtraction",
    round_to_int: bool = False,
    p_cut: float = 0.01,
    verbose: bool = True,
) -> AnnData:
    """Remove ambient RNA contamination via Rust core.

    Writes corrected counts to adata.layers["soupx_corrected"].
    Modifies adata in-place and returns it for chaining.

    Parameters
    ----------
    adata
        AnnData with cell counts and soupx_rho in adata.obs.
    method
        "subtraction" (default), "soup_only", or "multinomial".
    round_to_int
        If True, round corrected counts to integers.
    p_cut
        p-value cutoff for soup_only method.
    verbose
        If True, print progress messages.

    Returns
    -------
    The same AnnData object (in-place + chaining).
    """
    _ensure_core()

    if "soupx_rho" not in adata.obs:
        raise ValueError(
            "adata.obs['soupx_rho'] not found. Run auto_est_cont() or "
            "estimate_contamination() first."
        )

    rho = float(adata.obs["soupx_rho"].iloc[0] if hasattr(adata.obs["soupx_rho"], 'iloc') else adata.obs["soupx_rho"][0])

    if verbose:
        print(f"adjust_counts: method={method}, ρ={rho:.4f}")

    # Prepare CSR input
    toc = sparse.csr_matrix(adata.X) if not sparse.issparse(adata.X) else adata.X.tocsr()
    data, indices, indptr, shape = _csr_to_rust_args(toc)
    soup_profile = adata.uns["soupx_soup_profile"]

    # Call Rust
    res_data, res_indices, res_indptr, (n_rows, n_cols) = _core.adjust_counts(
        data, indices, indptr, shape,
        soup_profile.tolist(),
        rho,
        method,
        p_cut,
    )

    # Reconstruct scipy CSR
    n_rows = int(n_rows)
    n_cols = int(n_cols)
    corrected = sparse.csr_matrix(
        (np.array(res_data), np.array(res_indices), np.array(res_indptr)),
        shape=(n_rows, n_cols),
    )

    if round_to_int:
        corrected.data = np.round(corrected.data).astype(np.int32)
        corrected.eliminate_zeros()

    adata.layers["soupx_corrected"] = corrected
    adata.uns["soupx_method"] = method

    if verbose:
        nnz_orig = toc.nnz
        nnz_corr = corrected.nnz
        print(f"  Original: {nnz_orig} nnz → Corrected: {nnz_corr} nnz "
              f"({100 * (1 - nnz_corr / max(nnz_orig, 1)):.1f}% reduction)")

    return adata


# ============================================================================
# decontaminate — convenience wrapper
# ============================================================================

def decontaminate(
    adata: AnnData,
    *,
    clusters: str = "leiden",
    method: Literal["subtraction", "soup_only", "multinomial"] = "subtraction",
    round_to_int: bool = True,
    **auto_est_kwargs,
) -> AnnData:
    """Convenience: set clusters → auto-estimate ρ → adjust counts.

    Prerequisites: load_10x() already called, scanpy clustering done.
    Modifies adata in-place and returns it for chaining.
    """
    set_clusters(adata, key=clusters)
    auto_est_cont(adata, **auto_est_kwargs)
    adjust_counts(adata, method=method, round_to_int=round_to_int)
    return adata


# ============================================================================
# quick_markers — calls Rust core
# ============================================================================

def quick_markers(
    adata: AnnData,
    *,
    n: int = 10,
    fdr: float = 0.01,
    express_cut: float = 0.9,
) -> pd.DataFrame:
    """Quick marker gene discovery via tf-idf + hypergeometric test (Rust core).

    Parameters
    ----------
    adata
        AnnData with cell counts.
    n
        Number of top marker genes per cluster.
    fdr
        False discovery rate threshold.
    express_cut
        Binarization cutoff (counts > cut → expressed).

    Returns
    -------
    DataFrame with columns: cluster, gene, gene_idx, tfidf_score, qval.
    """
    _ensure_core()

    toc = sparse.csr_matrix(adata.X) if not sparse.issparse(adata.X) else adata.X.tocsr()
    data, indices, indptr, shape = _csr_to_rust_args(toc)
    clusters_vec = _clusters_to_rust(adata).tolist()
    var_names = list(adata.var_names)

    results = _core.quick_markers(
        data, indices, indptr, shape,
        clusters_vec, n, fdr, express_cut,
    )

    # Build DataFrame
    rows = []
    for cluster_id, cluster_markers in enumerate(results):
        for entry in cluster_markers:
            gene_idx, tfidf, qval = entry
            rows.append({
                "cluster": cluster_id,
                "gene_idx": gene_idx,
                "gene": var_names[gene_idx] if gene_idx < len(var_names) else str(gene_idx),
                "tfidf_score": tfidf,
                "qval": qval,
            })

    return pd.DataFrame(rows)


# ============================================================================
# 10X CellRanger file parsing (pure Python)
# ============================================================================

def _find_raw_matrix_dir(cellranger_dir: Path) -> Path:
    """Find raw_feature_bc_matrix directory in CellRanger output."""
    candidates = [
        cellranger_dir / "raw_feature_bc_matrix",
        cellranger_dir / "raw_gene_bc_matrices",
        cellranger_dir / "raw_gene_bc_matrix",
    ]
    for c in candidates:
        if _is_matrix_dir(c):
            return c
    for pattern in ["raw_feature_bc_matrix", "raw_gene_bc_matrices"]:
        for m in cellranger_dir.rglob(pattern):
            if _is_matrix_dir(m):
                return m
    raise FileNotFoundError(
        f"Could not find raw matrix in {cellranger_dir}"
    )


def _is_matrix_dir(path: Path) -> bool:
    return path.is_dir() and (
        (path / "matrix.mtx.gz").exists() or (path / "matrix.mtx").exists()
    )


def _read_mtx_to_csr(matrix_dir: Path) -> sparse.csr_matrix:
    """Read 10X mtx directory into CSR (genes × barcodes)."""
    from scipy.io import mmread

    mtx_path = _find_file(matrix_dir, "matrix.mtx")
    matrix = sparse.csr_matrix(
        mmread(gzip.open(mtx_path, "rt") if mtx_path.suffix == ".gz" else str(mtx_path)).tocoo()
    )

    barcodes = _read_tsv_col(_find_file(matrix_dir, "barcodes.tsv"))
    gene_ids, gene_names = _read_features(_find_file(matrix_dir, "features.tsv", "genes.tsv"))

    if matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(gene_ids):
        matrix = matrix.T.tocsr()

    return matrix.astype(np.float64)


def _find_file(directory: Path, *names: str) -> Path:
    for name in names:
        for variant in (name + ".gz", name):
            p = directory / variant
            if p.exists():
                return p
    raise FileNotFoundError(f"Could not find any of {names} in {directory}")


def _read_tsv_col(path: Path, col: int = 0) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.strip().split("\t")[col] for line in f if line.strip()]


def _read_features(path: Path) -> tuple[list[str], list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    ids, names = [], []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            ids.append(parts[0])
            names.append(parts[1] if len(parts) >= 2 else parts[0])
    return ids, names


def _store_soup_profile(
    adata: AnnData,
    soup_profile: np.ndarray,
    soup_counts: np.ndarray,
    soup_range: tuple[int, int] = (0, 100),
    n_empty: int = 0,
) -> None:
    """Store soup profile in AnnData."""
    adata.uns["soupx_soup_profile"] = soup_profile
    adata.varm["soup_profile_est"] = soup_profile.reshape(-1, 1)
    # soup_counts may be from full gene set; only store if aligned
    if len(soup_counts) == len(soup_profile):
        adata.varm["soup_profile_counts"] = soup_counts.reshape(-1, 1)
        adata.uns["soupx_total_soup_umis"] = float(soup_counts.sum())
    adata.uns["soupx_soup_range"] = list(soup_range)
    adata.uns["soupx_n_empty"] = n_empty


def _count_lines(path: Path) -> int:
    """Count lines in a file (supports .gz)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return sum(1 for _ in f)
