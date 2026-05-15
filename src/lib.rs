use pyo3::prelude::*;

pub mod sparse;
pub mod stats;
pub mod soup_profile;
pub mod markers;
pub mod non_expressing;
pub mod contamination;
pub mod adjustment;
pub mod alloc;
pub mod expand;

/// soupX Python module (Rust core engine).
///
/// Ambient RNA contamination removal for single-cell RNA-seq data,
/// reimplemented from the R SoupX package.
#[pymodule]
#[pyo3(name = "soupx_core")]
mod soupx {
    use super::*;
    use numpy::PyArray1;
    use pyo3::types::{PyList, PyTuple};
    use pyo3::Bound;

    // ========================================================================
    // Utility
    // ========================================================================

    /// Return the crate version.
    #[pyfunction]
    fn _version() -> PyResult<String> {
        Ok(env!("CARGO_PKG_VERSION").to_string())
    }

    /// Validate sparse matrix bridge: returns (nnz, n_rows, n_cols).
    #[pyfunction]
    fn _validate_csr(
        data: Bound<'_, PyArray1<f64>>,
        indices: Bound<'_, PyArray1<i64>>,
        indptr: Bound<'_, PyArray1<i64>>,
        shape: (usize, usize),
    ) -> PyResult<(usize, usize, usize)> {
        let csr = super::sparse::csr_from_numpy_owned(&data, &indices, &indptr, shape)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok((csr.nnz(), csr.rows(), csr.cols()))
    }

    // ========================================================================
    // quick_markers
    // ========================================================================

    /// Quick marker gene discovery via tf-idf + hypergeometric test.
    ///
    /// Returns a list of lists: one list per cluster, each entry is
    /// [gene_idx, tfidf_score, qval].
    #[pyfunction]
    fn quick_markers(
        py: Python<'_>,
        data: Bound<'_, PyArray1<f64>>,
        indices: Bound<'_, PyArray1<i64>>,
        indptr: Bound<'_, PyArray1<i64>>,
        shape: (usize, usize),
        clusters: Vec<usize>,
        n: usize,
        fdr: f64,
        express_cut: f64,
    ) -> PyResult<Py<PyList>> {
        let toc = super::sparse::csr_from_numpy_owned(&data, &indices, &indptr, shape)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let results = super::markers::quick_markers(&toc, &clusters, n, fdr, express_cut);

        let py_results = PyList::empty(py);
        for cluster_markers in &results {
            let py_cluster = PyList::empty(py);
            for m in cluster_markers {
                let entry = PyList::empty(py);
                entry.append(m.gene_idx)?;
                entry.append(m.tfidf_score)?;
                entry.append(m.qval)?;
                py_cluster.append(entry)?;
            }
            py_results.append(py_cluster)?;
        }

        Ok(py_results.unbind())
    }

    // ========================================================================
    // auto_est_cont
    // ========================================================================

    /// Automatically estimate contamination fraction ρ.
    ///
    /// Returns estimated ρ as a float, or raises ValueError if estimation fails.
    #[pyfunction]
    fn auto_est_cont(
        data: Bound<'_, PyArray1<f64>>,
        indices: Bound<'_, PyArray1<i64>>,
        indptr: Bound<'_, PyArray1<i64>>,
        shape: (usize, usize),
        clusters: Vec<usize>,
        soup_profile: Vec<f64>,
        tfidf_min: f64,
        soup_quantile: f64,
        max_markers: usize,
        contamination_range: (f64, f64),
        prior_rho: f64,
        prior_rho_stddev: f64,
    ) -> PyResult<f64> {
        let toc = super::sparse::csr_from_numpy_owned(&data, &indices, &indptr, shape)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let rho = super::contamination::auto_est_cont(
            &toc,
            &clusters,
            &soup_profile,
            tfidf_min,
            soup_quantile,
            max_markers,
            contamination_range,
            prior_rho,
            prior_rho_stddev,
        );

        match rho {
            Some(r) => Ok(r),
            None => Err(pyo3::exceptions::PyValueError::new_err(
                "auto_est_cont: estimation failed — no valid probes found"
            )),
        }
    }

    // ========================================================================
    // estimate_contamination (manual mode)
    // ========================================================================

    /// Manual contamination estimation from user-specified gene sets.
    ///
    /// `gene_sets` is a list of lists of global gene indices.
    #[pyfunction]
    fn estimate_contamination(
        data: Bound<'_, PyArray1<f64>>,
        indices: Bound<'_, PyArray1<i64>>,
        indptr: Bound<'_, PyArray1<i64>>,
        shape: (usize, usize),
        gene_sets: Vec<Vec<usize>>,
        clusters: Vec<usize>,
        soup_profile: Vec<f64>,
        max_contam: f64,
        fdr: f64,
    ) -> PyResult<f64> {
        let toc = super::sparse::csr_from_numpy_owned(&data, &indices, &indptr, shape)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let rho = super::contamination::calculate_contamination_fraction(
            &toc,
            &gene_sets,
            &clusters,
            &soup_profile,
            max_contam,
            fdr,
        );

        match rho {
            Some(r) => Ok(r),
            None => Err(pyo3::exceptions::PyValueError::new_err(
                "estimate_contamination: estimation failed"
            )),
        }
    }

    // ========================================================================
    // adjust_counts
    // ========================================================================

    /// Remove ambient RNA contamination from cell counts.
    ///
    /// Returns a tuple (data, indices, indptr, (n_rows, n_cols)) representing
    /// the corrected CSR matrix. Reconstruct in Python with:
    ///   scipy.sparse.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
    #[pyfunction]
    fn adjust_counts(
        py: Python<'_>,
        data: Bound<'_, PyArray1<f64>>,
        indices: Bound<'_, PyArray1<i64>>,
        indptr: Bound<'_, PyArray1<i64>>,
        shape: (usize, usize),
        soup_profile: Vec<f64>,
        rho: f64,
        method: String,
        p_cut: f64,
    ) -> PyResult<Py<PyTuple>> {
        let toc = super::sparse::csr_from_numpy_owned(&data, &indices, &indptr, shape)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let result = super::adjustment::adjust_counts(
            &toc,
            &soup_profile,
            rho,
            &method,
            p_cut,
        );

        // Build numpy arrays for data, indices, indptr
        let n_rows = result.rows();
        let n_cols = result.cols();
        let (result_indptr_raw, result_indices_raw, result_data) = result.into_raw_storage();

        let py_data = PyArray1::from_vec(py, result_data);
        // Convert usize indices to i64 for scipy compatibility
        let indices_i64: Vec<i64> = result_indices_raw.iter().map(|&i| i as i64).collect();
        let indptr_i64: Vec<i64> = result_indptr_raw.iter().map(|&i| i as i64).collect();
        let py_indices = PyArray1::from_vec(py, indices_i64);
        let py_indptr = PyArray1::from_vec(py, indptr_i64);

        // Build result tuple — each element as Py<PyAny>
        let data_obj: Py<PyAny> = py_data.into_pyobject(py)?.into_any().unbind();
        let indices_obj: Py<PyAny> = py_indices.into_pyobject(py)?.into_any().unbind();
        let indptr_obj: Py<PyAny> = py_indptr.into_pyobject(py)?.into_any().unbind();
        let shape_obj: Py<PyAny> = PyTuple::new(py, [
            (n_rows as i64).into_pyobject(py)?,
            (n_cols as i64).into_pyobject(py)?,
        ])?.into_any().unbind();

        Ok(PyTuple::new(py, [data_obj, indices_obj, indptr_obj, shape_obj])?.unbind())
    }

    // ========================================================================
    // BH correction (utility)
    // ========================================================================

    /// Benjamini-Hochberg FDR correction.
    /// Returns q-values (same length as input).
    #[pyfunction]
    fn bh_correct(p_values: Vec<f64>) -> Vec<f64> {
        super::stats::bh_correction(&p_values)
    }
}
