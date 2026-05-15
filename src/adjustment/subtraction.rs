//! Subtraction method for count adjustment.
//!
//! Iteratively subtracts expected background counts from each gene per cell,
//! using the alloc algorithm to distribute the total soup UMI budget across genes.
//!
//! Uses COO triplets as intermediate format to correctly build CSR output
//! regardless of per-cell processing order.
//!
//! See REFACTOR_PLAN.md §1.4f for algorithm details.

use crate::alloc;

/// Subtract background contamination from cell counts.
///
/// # Arguments
/// * `toc` - Cell count matrix (CSR, genes × cells)
/// * `soup_profile` - Estimated soup proportions per gene (global, length = n_genes)
/// * `rho` - Global contamination fraction
///
/// # Returns
/// Corrected count matrix in CSR format (genes × cells).
pub fn subtract(
    toc: &sprs::CsMatI<f64, usize>,
    soup_profile: &[f64],
    rho: f64,
) -> sprs::CsMatI<f64, usize> {
    let (n_genes, n_cells) = (toc.rows(), toc.cols());
    let tol = 1e-6;
    let max_iter = 100;

    // Build result as COO triplets (row, col, value), then convert to CSR.
    let mut coo_rows: Vec<usize> = Vec::new();
    let mut coo_cols: Vec<usize> = Vec::new();
    let mut coo_data: Vec<f64> = Vec::new();

    for cell_idx in 0..n_cells {
        let (col_data, gene_indices) = get_column(toc, cell_idx);
        let n_umis: f64 = col_data.iter().sum();
        let exp_soup = n_umis * rho;

        let mut fit: Vec<f64> = col_data.clone();
        let mut remaining = exp_soup;
        let mut iter = 0;

        while remaining > tol && iter < max_iter {
            let allocated = alloc::alloc(
                remaining,
                &fit,
                soup_profile,
                &gene_indices,
            );
            for (i, &a) in allocated.iter().enumerate() {
                fit[i] = (fit[i] - a).max(0.0);
            }
            remaining -= allocated.iter().sum::<f64>();
            iter += 1;
        }

        // Emit COO triplets for non-zero corrected counts
        for (i, &val) in fit.iter().enumerate() {
            if val > 0.0 {
                coo_rows.push(gene_indices[i]); // global gene index = row
                coo_cols.push(cell_idx);
                coo_data.push(val);
            }
        }
    }

    // Build CSR from COO triplets
    coo_to_csr(n_genes, n_cells, &coo_rows, &coo_cols, &coo_data)
}

/// Convert COO (triplet) format to CSR.
///
/// COO triplets are sorted by row then column; duplicate entries are summed.
fn coo_to_csr(
    n_rows: usize,
    n_cols: usize,
    coo_rows: &[usize],
    coo_cols: &[usize],
    coo_data: &[f64],
) -> sprs::CsMatI<f64, usize> {
    let nnz = coo_data.len();

    // Sort by (row, col)
    let mut indexed: Vec<(usize, usize, f64)> = (0..nnz)
        .map(|i| (coo_rows[i], coo_cols[i], coo_data[i]))
        .collect();
    indexed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    // Build CSR: indptr has n_rows + 1 elements
    let mut indptr = vec![0usize; n_rows + 1];
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);

    let mut current_row = 0;
    for (row, col, val) in &indexed {
        while current_row < *row {
            current_row += 1;
            indptr[current_row] = data.len();
        }
        // Merge duplicates (same row & col) by summing
        if !data.is_empty()
            && indices.last() == Some(col)
            && current_row == *row
        {
            *data.last_mut().unwrap() += val;
        } else {
            indices.push(*col);
            data.push(*val);
        }
    }
    // Fill remaining indptr entries
    for r in (current_row + 1)..=n_rows {
        indptr[r] = data.len();
    }

    sprs::CsMatI::new((n_rows, n_cols), indptr, indices, data)
}

/// Extract a single column from a CSR matrix.
///
/// Returns (values, global_gene_indices) for non-zero entries in the column.
///
/// Complexity: O(nnz) scan per cell. For production, pre-build CSC or use
/// indptr chunking + column index. See REFACTOR_PLAN.md §4.4.
fn get_column(mat: &sprs::CsMatI<f64, usize>, col: usize) -> (Vec<f64>, Vec<usize>) {
    let mut data = Vec::new();
    let mut indices = Vec::new();
    for row in 0..mat.rows() {
        if let Some(val) = mat.get(row, col) {
            data.push(*val);
            indices.push(row);
        }
    }
    (data, indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coo_to_csr_empty() {
        let csr = coo_to_csr(2, 3, &[], &[], &[]);
        assert_eq!(csr.rows(), 2);
        assert_eq!(csr.cols(), 3);
        assert_eq!(csr.nnz(), 0);
    }

    #[test]
    fn test_coo_to_csr_basic() {
        let csr = coo_to_csr(
            3,
            2,
            &[0, 0, 1],
            &[0, 1, 0],
            &[1.0, 2.0, 3.0],
        );
        assert_eq!(csr.rows(), 3);
        assert_eq!(csr.cols(), 2);
        assert_eq!(csr.nnz(), 3);
        assert_eq!(*csr.get(0, 0).unwrap(), 1.0);
        assert_eq!(*csr.get(0, 1).unwrap(), 2.0);
        assert_eq!(*csr.get(1, 0).unwrap(), 3.0);
    }

    #[test]
    fn test_coo_to_csr_merge_duplicates() {
        // Two entries for (0, 0): 1.0 + 2.0 = 3.0
        let csr = coo_to_csr(
            2,
            2,
            &[0, 0],
            &[0, 0],
            &[1.0, 2.0],
        );
        assert_eq!(csr.nnz(), 1);
        assert_eq!(*csr.get(0, 0).unwrap(), 3.0);
    }

    #[test]
    fn test_subtract_zero_rho() {
        // With rho=0, no subtraction → output = input
        let indptr = vec![0usize, 2, 3, 4];
        let indices = vec![0usize, 1, 0, 1];
        let data = vec![5.0, 3.0, 2.0, 1.0];
        let toc = sprs::CsMatI::new((3, 2), indptr, indices, data);
        let profile = vec![1.0 / 3.0; 3];

        let result = subtract(&toc, &profile, 0.0);
        assert_eq!(result.rows(), 3);
        assert_eq!(result.cols(), 2);
        // With rho=0, all entries should be preserved
        assert_eq!(result.nnz(), 4);
    }
}
