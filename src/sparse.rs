//! CSR/CSC sparse matrix bridge between Python numpy/scipy and Rust sprs.
//!
//! Strategy (safety-first):
//! - Default: clone data from numpy buffers into owned `sprs::CsMat` in GIL-held context,
//!   then release GIL for all subsequent computation.
//! - Alternative (unsafe, opt-in): zero-copy CsMatView for short-lived single-call operations
//!   that never cross `py.allow_threads()` boundaries.
//!
//! See REFACTOR_PLAN.md §4.3(a) for rationale.

use numpy::PyArray1;
use numpy::PyArrayMethods as _;
use pyo3::Bound;

/// Clone a scipy CSR matrix from numpy arrays into an owned sprs::CsMat.
///
/// This is the **safe default** path: data is copied during GIL-held context,
/// after which Rust fully owns the matrix and can freely use rayon + allow_threads.
///
/// Typical clone cost: 50K cells × 25K genes × 5% sparsity ≈ 0.15 GB, ~0.1s.
pub fn csr_from_numpy_owned(
    data: &Bound<'_, PyArray1<f64>>,
    indices: &Bound<'_, PyArray1<i64>>,
    indptr: &Bound<'_, PyArray1<i64>>,
    shape: (usize, usize),
) -> Result<sprs::CsMatI<f64, usize>, Error> {
    let data_readonly = data.readonly();
    let data_slice = data_readonly
        .as_slice()
        .map_err(|_| Error::NumpyAccessError)?;

    let indices_readonly = indices.readonly();
    let indices_slice = indices_readonly
        .as_slice()
        .map_err(|_| Error::NumpyAccessError)?;

    let indptr_readonly = indptr.readonly();
    let indptr_slice = indptr_readonly
        .as_slice()
        .map_err(|_| Error::NumpyAccessError)?;

    let indices_usize: Vec<usize> = indices_slice.iter().map(|&i| i as usize).collect();
    let indptr_usize: Vec<usize> = indptr_slice.iter().map(|&i| i as usize).collect();

    Ok(sprs::CsMatI::new(
        shape,
        indptr_usize,
        indices_usize,
        data_slice.to_vec(),
    ))
}

/// Error type for sparse matrix operations.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Failed to access numpy array data")]
    NumpyAccessError,
    #[error("Invalid sparse matrix shape or dimensions")]
    InvalidShape,
}
