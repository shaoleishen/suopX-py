"""Tests for soupX Python + Rust package."""

# TODO: comprehensive tests
# - Level 0: statistical function accuracy vs R
# - Level 1: soup profile estimation
# - Level 2: quickMarkers vs R reference
# - Level 3: ρ estimation accuracy
# - Level 4: corrected count matrix consistency
# See REFACTOR_PLAN.md §8 "对账策略"


def test_import():
    """Basic import test."""
    import soupx
    assert soupx.__version__ == "0.1.0"


def test_soup_channel():
    """Test SoupChannel construction."""
    import numpy as np
    import scipy.sparse as sparse
    from soupx import SoupChannel

    toc = sparse.csr_matrix(np.eye(3))
    tod = sparse.csr_matrix(np.ones((3, 10)))
    # With explicit soup_profile to avoid estimate_soup
    sc = SoupChannel(toc, tod, soup_profile=np.ones(3) / 3.0)
    assert sc.n_genes == 3
    assert sc.n_cells == 3
    assert sc.n_droplets == 10


def test_set_clusters():
    """Test set_clusters."""
    import numpy as np
    import anndata
    from soupx import set_clusters

    adata = anndata.AnnData(np.eye(5))
    adata.obs["test_clusters"] = ["A", "A", "B", "B", "C"]
    set_clusters(adata, key="test_clusters")
    assert adata.uns["soupx_cluster_key"] == "test_clusters"


def test_decontaminate_requires_load():
    """Test decontaminate raises if soup data not loaded."""
    import numpy as np
    import anndata
    from soupx import decontaminate

    adata = anndata.AnnData(np.eye(5))
    adata.obs["leiden"] = ["0", "0", "1", "1", "2"]
    adata.uns["soupx_cluster_key"] = "leiden"
    # Missing soupx_soup_profile → KeyError
    import pytest
    with pytest.raises((KeyError, ImportError)):
        decontaminate(adata)
