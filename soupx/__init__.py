# soupx
#
# Python package for ambient RNA contamination removal in single-cell RNA-seq data.
# Rust core engine, Python API integrated with scanpy/AnnData ecosystem.
#
# Original R package: constantAmateur/SoupX (GigaScience 2020)
# Python + Rust reimplementation: suopX-py

# ---------------------
# Version
# ---------------------
from soupx._version import __version__

# ---------------------
# Public API
# ---------------------
from soupx._bridge import (
    load_10x,
    set_clusters,
    auto_est_cont,
    estimate_contamination,
    adjust_counts,
    decontaminate,
    quick_markers,
)
from soupx._channel import SoupChannel

# ---------------------
# Visualization (import on demand)
# ---------------------
from soupx import plotting

__all__ = [
    "__version__",
    # Core API
    "load_10x",
    "set_clusters",
    "auto_est_cont",
    "estimate_contamination",
    "adjust_counts",
    "decontaminate",
    "quick_markers",
    # Low-level
    "SoupChannel",
    # Submodules
    "plotting",
]
