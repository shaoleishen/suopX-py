//! Expand cluster-level results to per-cell results.
//!
//! Some algorithms (e.g., adjustCounts in R) operate on cluster-aggregated
//! matrices and then expand back to per-cell matrices.
//!
//! Parallelized across clusters with rayon.

/// Expand cluster-level matrix to per-cell matrix.
///
/// # Arguments
/// * `cluster_mat` - Result per cluster (genes × n_clusters)
/// * `clusters` - Cluster assignment per cell (length = n_cells)
///
/// # Returns
/// Per-cell matrix (genes × n_cells)
pub fn expand_clusters(
    _cluster_mat: &sprs::CsMatI<f64, usize>,
    _clusters: &[usize],
) -> sprs::CsMatI<f64, usize> {
    // Placeholder - TODO: implement
    unimplemented!("cluster-to-cell expansion")
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_empty() {
        // TODO
    }
}
