#!/usr/bin/env python3
"""Benchmark: soupX-py on real 10X snRNA-seq data (M7 sample).

468K raw droplets, 16K cells, 39K genes.
"""

import os, sys, time, resource, gc
import numpy as np
import scipy.sparse as sparse
import scanpy as sc

SAMPLE = "M7"
DATA_DIR = f"/home/bioshen/CLI/snRNAseq/{SAMPLE}"

def mem_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)

def log(msg, start=None):
    e = f" [+{time.time()-start:.1f}s]" if start else ""
    print(f"[{mem_gb():.2f}GB] {msg}{e}")

t_total = time.time()

# ---------------------------------------------------------------
# 1. Load filtered data
# ---------------------------------------------------------------
t0 = time.time()
log("Loading filtered matrix...", t0)
adata = sc.read_10x_mtx(DATA_DIR + "/filtered_feature_bc_matrix/", gex_only=False)
log(f"  Shape: {adata.shape} | nnz: {adata.X.nnz}", t0)

# Quick clustering for demo
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=15)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.5)
log(f"  Clusters: {adata.obs['leiden'].nunique()}", t0)

# ---------------------------------------------------------------
# 2. soupX-py: load_10x (soup profile from raw droplets)
# ---------------------------------------------------------------
t1 = time.time()
log("=== load_10x ===", t1)
import soupx
soupx.load_10x(adata, cellranger_dir=DATA_DIR, soup_range=(0, 100), keep_droplets=False)
n_empty = adata.uns.get("soupx_n_empty", 0)
log(f"  Empty droplets used: {n_empty}", t1)

# ---------------------------------------------------------------
# 3. soupX-py: auto_est_cont
# ---------------------------------------------------------------
t2 = time.time()
log("=== auto_est_cont ===", t2)
soupx.set_clusters(adata, key="leiden")
soupx.auto_est_cont(adata, tfidf_min=1.0, max_markers=100, verbose=False)
rho = adata.obs["soupx_rho"].iloc[0]
log(f"  Estimated ρ = {rho:.4f}", t2)

# ---------------------------------------------------------------
# 4. soupX-py: adjust_counts (subtraction only)
# ---------------------------------------------------------------
t3 = time.time()
log("=== adjust_counts (subtraction) ===", t3)
soupx.adjust_counts(adata, method="subtraction", verbose=False)
corr_nnz = adata.layers["soupx_corrected"].nnz
orig_nnz = adata.X.nnz
reduction = 100 * (1 - corr_nnz / max(orig_nnz, 1))
log(f"  nnz: {orig_nnz} → {corr_nnz} ({reduction:.1f}% reduction)", t3)

# ---------------------------------------------------------------
# 5. soupX-py: adjust_counts (soup_only)
# ---------------------------------------------------------------
t4 = time.time()
log("=== adjust_counts (soup_only) ===", t4)
soupx.adjust_counts(adata, method="soup_only", p_cut=0.01, verbose=False)
log(f"  Done", t4)

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
gc.collect()
print()
print("=" * 60)
print("  BENCHMARK: soupX-py on M7")
print("=" * 60)
print(f"  Cells: {adata.shape[0]} | Genes: {adata.shape[1]}")
print(f"  Raw droplets: {adata.uns.get('soupx_n_droplets', '?')}")
print(f"  Empty droplets: {n_empty}")
print(f"  Clusters: {adata.obs['leiden'].nunique()}")
print(f"  ρ (auto): {rho:.4f}")
print(f"  Total wall time: {time.time() - t_total:.1f}s")
print(f"  Peak RSS: {mem_gb():.2f} GB")
print("=" * 60)
