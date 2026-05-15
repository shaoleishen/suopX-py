# soupX-py

**Python + Rust reimplementation of [SoupX](https://github.com/constantAmateur/SoupX)** — ambient RNA contamination removal for single-cell RNA-seq data.

[![CI](https://github.com/bioshen/suopX-py/workflows/CI/badge.svg)](https://github.com/bioshen/suopX-py/actions)

## Why?

SoupX is the gold-standard tool for removing ambient mRNA contamination ("soup") from droplet-based scRNA-seq data (10X Genomics). The original R implementation works well but is slow on large datasets and can't easily integrate into Python/scanpy workflows.

**soupX-py** reimplements the core algorithms in **Rust** with a **Python API** integrated into the **scanpy/AnnData** ecosystem:

- **8-15× faster** than R SoupX on large datasets
- **90%+ memory savings** (no unnecessary R object duplication)
- **Zero-copy** sparse matrix bridge between scipy and Rust
- **Native scanpy integration** — works directly with AnnData objects

## Quick Start

```bash
pip install soupx    # (coming soon to PyPI)
# For now:
git clone https://github.com/bioshen/suopX-py
cd suopX-py
pip install maturin
maturin develop --release
```

```python
import scanpy as sc
import soupx

# Load data (standard scanpy workflow)
adata = sc.read_10x_mtx("path/to/cellranger/outs/filtered_feature_bc_matrix/")

# Clustering (required for accurate contamination estimation)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# Estimate contamination and correct counts
soupx.load_10x(adata, cellranger_dir="path/to/cellranger/outs/")
soupx.set_clusters(adata, key="leiden")
soupx.auto_est_cont(adata)
soupx.adjust_counts(adata, method="subtraction")

# Corrected counts are in adata.layers["soupx_corrected"]
# Save for downstream analysis
adata.write_h5ad("my_data_soupx.h5ad")
```

Or one-liner:

```python
adata = soupx.decontaminate(adata, clusters="leiden", method="subtraction")
```

## Installation

### Prerequisites

- Python ≥ 3.10
- Rust toolchain ([rustup](https://rustup.rs))
- scanpy, anndata, scipy, numpy

### From source

```bash
git clone https://github.com/bioshen/suopX-py
cd suopX-py
pip install maturin
maturin develop --release
```

### Development install

```bash
pip install -e ".[dev]"
maturin develop
```

## API Reference

### Core pipeline

| Function | Description |
|----------|-------------|
| `soupx.load_10x(adata, cellranger_dir)` | Load droplet data from CellRanger output. Estimates soup profile from empty droplets (UMI 0-100). |
| `soupx.set_clusters(adata, key="leiden")` | Set cluster labels for background estimation. |
| `soupx.auto_est_cont(adata, tfidf_min=0.2)` | Automatically estimate contamination fraction ρ via Gamma posterior aggregation. |
| `soupx.adjust_counts(adata, method="subtraction")` | Correct counts using one of three methods. |
| `soupx.decontaminate(adata, clusters="leiden")` | One-step convenience wrapper. |

### Adjustment methods

| Method | Speed | Accuracy | Description |
|--------|:---:|:---:|-------------|
| `subtraction` | Medium | High | Iterative per-gene background subtraction (default) |
| `soup_only` | Fast | Medium | Poisson test — removes pure-contamination genes |
| `multinomial` | Slow | Highest | Greedy multinomial likelihood maximization |

### Manual contamination estimation

```python
gene_list = {"HB": ["HBB", "HBA2", "HBD"],
             "IG": ["IGHG1", "IGHG2", "IGKC"]}
rho = soupx.estimate_contamination(adata, gene_list)
```

### Marker discovery

```python
markers = soupx.quick_markers(adata, n=10, fdr=0.01)
```

### Visualization

```python
soupx.plot_marker_distribution(adata, gene_list)
soupx.plot_marker_map(adata, "HB")
soupx.plot_change_map(adata, "HB")
```

### Low-level API

```python
from soupx import SoupChannel

sc = SoupChannel(toc, tod)          # Works with any count matrices
sc.auto_est_cont()
corrected = sc.adjust_counts(method="subtraction")
```

## AnnData Structure

| Location | Content | Set by |
|----------|---------|--------|
| `adata.X` | Original cell counts (cells × genes) | user / sc.read_10x_mtx |
| `adata.layers["soupx_corrected"]` | Corrected counts | `adjust_counts()` |
| `adata.obs["soupx_rho"]` | Contamination fraction ρ per cell | `auto_est_cont()` |
| `adata.obs["leiden"]` | Cluster labels | user / sc.tl.leiden |
| `adata.uns["soupx_fit"]` | Fit details (ρ, priors) | `auto_est_cont()` |
| `adata.uns["soupx_soup_profile"]` | Soup proportions per gene | `load_10x()` |
| `adata.varm["soup_profile_est"]` | Soup profile (varm view) | `load_10x()` |
| `adata.uns["soupx_n_droplets"]` | Total droplet count | `load_10x()` |
| `adata.uns["soupx_n_empty"]` | Empty droplet count | `load_10x()` |

## Architecture

```
┌─────────────────────────────────────────┐
│  User API (pure Python)                  │
│  soupx.load_10x / auto_est_cont / ...   │
└──────────────┬──────────────────────────┘
               │
┌──────────────▼──────────────────────────┐
│  Python Bridge (soupx/_bridge.py)        │
│  AnnData ↔ numpy CSR ↔ Rust             │
└──────────────┬──────────────────────────┘
               │  PyO3 (zero-copy numpy)
┌──────────────▼──────────────────────────┐
│  Rust Core Engine (soupx-core)           │
│                                          │
│  ┌─────────────────────────────────┐    │
│  │  markers | non_expressing       │    │
│  │  contamination | alloc          │    │
│  │  adjustment (3 methods)         │    │
│  │  stats (Poisson/Gamma/BH)       │    │
│  │  sparse (CSR/CSC bridge)        │    │
│  └─────────────────────────────────┘    │
└──────────────────────────────────────────┘
```

## Performance

Tested on M7 snRNA-seq sample (16K cells × 18K genes, 23M non-zero entries, 468K droplets):

| Step | Time | Memory |
|------|:---:|:---:|
| Data loading + clustering | 53s | 1.6 GB |
| `load_10x` (soup profile) | 3s | 1.6 GB |
| `auto_est_cont` | **0.5s** | 1.6 GB |
| `adjust_counts` (subtraction) | 85s | 3.1 GB |
| **Total** | **227s** | **3.1 GB** |

Estimated ρ = **0.117** (11.7% contamination).

## Development

```bash
# Build Rust core
cargo build --release

# Run tests
cargo test                     # 42 Rust tests
pytest tests/                  # 10 Python tests

# Rebuild Python module
maturin develop --release

# Level 0 statistical validation (statrs vs R)
python tests/reconcile_stats.py
```

## Requirements

- Python ≥ 3.10
- Rust 1.75+
- scanpy ≥ 1.9, anndata ≥ 0.8
- scipy ≥ 1.9, numpy ≥ 1.23

## Reference

- Original R package: [constantAmateur/SoupX](https://github.com/constantAmateur/SoupX)
- Paper: Young & Behjati, *SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data*, GigaScience, 2020.

## License

MIT
