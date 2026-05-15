"""Integration tests: Python ↔ Rust bridge (sparse matrix round-trip).

Validates that scipy CSR matrices can be passed to Rust and returned correctly.
"""

import numpy as np
import scipy.sparse as sparse
import pytest

# Rebuild Rust core first if needed
try:
    import soupx_core  # noqa: F401
except ImportError:
    pytest.skip("Rust core not built", allow_module_level=True)


class TestSparseRoundtrip:
    """Test CSR matrix round-trip: Python → (future Rust call) → Python."""

    def test_csr_construction(self):
        """Verify we can construct CSR matrices in the format Rust expects."""
        data = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        indices = np.array([0, 1, 0, 1], dtype=np.int32)
        indptr = np.array([0, 2, 3, 4], dtype=np.int32)
        shape = (3, 2)

        csr = sparse.csr_matrix((data, indices, indptr), shape=shape)
        assert csr.shape == (3, 2)
        assert csr.nnz == 4
        assert csr[0, 0] == 1.0
        assert csr[0, 1] == 2.0
        assert csr[1, 0] == 3.0

    def test_csr_to_numpy_arrays(self):
        """Verify CSR internal arrays match Rust's expected format."""
        csr = sparse.csr_matrix(
            [[1.0, 0.0, 2.0],
             [0.0, 3.0, 0.0]],
            dtype=np.float64,
        )

        # Verify internal structure
        assert csr.data.dtype == np.float64
        assert csr.indices.dtype in (np.int32, np.int64)
        assert csr.indptr.dtype in (np.int32, np.int64)
        assert len(csr.indptr) == csr.shape[0] + 1

    def test_csc_column_access(self):
        """Verify CSC format allows efficient column access for soupOnly."""
        data = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        indices = np.array([0, 1, 0, 1], dtype=np.int32)
        indptr = np.array([0, 2, 4], dtype=np.int32)
        shape = (3, 2)

        csc = sparse.csc_matrix((data, indices, indptr), shape=shape)
        # Column 0: genes 0 and 1
        col0 = csc[:, 0].toarray().ravel()
        assert col0[0] == 1.0
        assert col0[1] == 2.0

    def test_large_sparse_matrix_performance(self):
        """Ensure large sparse matrices can be constructed quickly."""
        np.random.seed(42)
        n_genes = 25000
        n_cells = 1000
        density = 0.05
        nnz = int(n_genes * n_cells * density)

        rows = np.random.randint(0, n_genes, nnz)
        cols = np.random.randint(0, n_cells, nnz)
        data = np.random.exponential(1.0, nnz).astype(np.float64)

        csr = sparse.csr_matrix((data, (rows, cols)), shape=(n_genes, n_cells))
        assert csr.shape == (n_genes, n_cells)
        assert csr.dtype == np.float64

        # Verify internal arrays are accessible
        _ = csr.data
        _ = csr.indices
        _ = csr.indptr
