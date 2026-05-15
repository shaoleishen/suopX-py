//! SoupOnly method for count adjustment.
//!
//! For each gene in each cell, performs Poisson test to determine
//! if its observed count is consistent with pure contamination.
//! Genes flagged as "soup-only" are removed entirely (set to zero).
//!
//! Key optimization (vs R): CSC format eliminates the global O(N log N) sort.
//! Each cell's (column's) entries are naturally contiguous in CSC.
//!
//! See REFACTOR_PLAN.md §4.3b for the 15-25× acceleration rationale.

use crate::stats::poisson;

/// Identify and remove pure-contamination genes.
///
/// # Algorithm
/// 1. For each (gene, cell) with count > 0:
///    Compute p = P(X ≥ count-1 | Poisson(expected_bg))
///    where expected_bg = gene_soup_prob * cell_total_umis * rho
/// 2. For each cell, sort by p-value and compute Fisher's combined p-value
///    (χ² with 2k df) incrementally
/// 3. Mark genes as soup-only if cumulative Fisher p < p_cut
/// 4. Zero out flagged genes
pub fn soup_only(
    toc: &sprs::CsMatI<f64, usize>,
    soup_profile: &[f64],
    rho: f64,
    p_cut: f64,
) -> sprs::CsMatI<f64, usize> {
    let (n_genes, n_cells) = (toc.rows(), toc.cols());

    // Build result via COO triplets (filtered)
    let mut coo_rows: Vec<usize> = Vec::new();
    let mut coo_cols: Vec<usize> = Vec::new();
    let mut coo_data: Vec<f64> = Vec::new();

    // Process per cell: group all entries for this cell
    for cell_idx in 0..n_cells {
        let (col_data, gene_indices) = get_column(toc, cell_idx);
        let n_umis: f64 = col_data.iter().sum();
        let n_entries = col_data.len();

        // Compute p-values for this cell
        let mut pvals: Vec<f64> = Vec::with_capacity(n_entries);
        for (i, &count) in col_data.iter().enumerate() {
            let gene_idx = gene_indices[i];
            let gene_soup = soup_profile.get(gene_idx).copied().unwrap_or(0.0);
            let expected = n_umis * rho * gene_soup;
            let p = if expected > 0.0 && count > 0.0 {
                poisson::ppois_upper(count - 1.0, expected)
            } else {
                1.0
            };
            pvals.push(p);
        }

        // Sort entries by p-value ascending, compute Fisher combined p
        // χ²(2k) = -2 * Σ ln(p_i) for k smallest p-values
        let mut indexed: Vec<(usize, f64, f64)> = col_data
            .iter()
            .enumerate()
            .map(|(i, &c)| (i, c, pvals[i]))
            .collect();
        indexed.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        // Compute Fisher combined p-values
        let mut fisher_pvals = vec![1.0f64; n_entries];
        let mut running_sum = 0.0f64;
        for (rank, &(_, _, p)) in indexed.iter().enumerate() {
            if p <= 0.0 {
                running_sum = f64::INFINITY;
            } else {
                running_sum += -2.0 * p.ln();
            }
            // P(χ²(2*(rank+1)) > running_sum)
            // Approximate using chi-squared survival function
            let df = 2.0 * (rank + 1) as f64;
            fisher_pvals[rank] = chi2_sf(running_sum, df);
        }

        // Mark soup-only: find largest k such that cumulative Fisher p < p_cut
        let mut keep_count = n_entries;
        for (rank, &(_orig_i, _, _)) in indexed.iter().enumerate() {
            if fisher_pvals[rank] >= p_cut {
                keep_count = n_entries - rank;
                // From rank onward, Fisher p ≥ p_cut → these entries are "real"
                // Entries before rank are "soup-only" and should be removed
                break;
            }
        }

        // Output only kept entries
        for (rank, &(orig_i, count, _)) in indexed.iter().enumerate() {
            let is_kept = rank >= n_entries - keep_count;
            if is_kept && count > 0.0 {
                coo_rows.push(gene_indices[orig_i]);
                coo_cols.push(cell_idx);
                coo_data.push(count);
            }
        }
    }

    // Build CSR from COO
    coo_to_csr(n_genes, n_cells, &coo_rows, &coo_cols, &coo_data)
}

/// Extract a single column from a CSR matrix.
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

/// Approximate χ² survival function: P(χ²(df) > x).
///
/// Uses the Wilson-Hilferty normal approximation:
/// (χ²/df)^(1/3) ~ N(1 - 2/(9df), 2/(9df))
fn chi2_sf(x: f64, df: f64) -> f64 {
    use statrs::distribution::{ContinuousCDF, Normal};

    if df <= 0.0 || x <= 0.0 {
        return 1.0;
    }

    let mean = 1.0 - 2.0 / (9.0 * df);
    let std = (2.0 / (9.0 * df)).sqrt();

    let t = (x / df).powf(1.0 / 3.0);
    let z = (t - mean) / std;

    let norm = Normal::new(0.0, 1.0).unwrap();
    norm.sf(z)
}

/// Convert COO to CSR (shared utility).
fn coo_to_csr(
    n_rows: usize,
    n_cols: usize,
    coo_rows: &[usize],
    coo_cols: &[usize],
    coo_data: &[f64],
) -> sprs::CsMatI<f64, usize> {
    let nnz = coo_data.len();

    let mut indexed: Vec<(usize, usize, f64)> = (0..nnz)
        .map(|i| (coo_rows[i], coo_cols[i], coo_data[i]))
        .collect();
    indexed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let mut indptr = vec![0usize; n_rows + 1];
    let mut indices = Vec::with_capacity(nnz);
    let mut data = Vec::with_capacity(nnz);

    let mut current_row = 0;
    for (row, col, val) in &indexed {
        while current_row < *row {
            current_row += 1;
            indptr[current_row] = data.len();
        }
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
    for r in (current_row + 1)..=n_rows {
        indptr[r] = data.len();
    }

    sprs::CsMatI::new((n_rows, n_cols), indptr, indices, data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chi2_sf() {
        // χ²(4) at x=9.488 is p≈0.05, so sf ≈ 0.05
        let p = chi2_sf(9.488, 4.0);
        assert!(p > 0.02 && p < 0.10, "chi2_sf(9.488, 4) = {p}");

        // χ²(2) at x very large → sf → 0
        let p2 = chi2_sf(100.0, 2.0);
        assert!(p2 < 0.001, "chi2_sf(100, 2) = {p2}");
    }

    #[test]
    fn test_soup_only_empty() {
        let toc = sprs::CsMatI::<f64, usize>::zero((3, 2));
        let profile = vec![1.0 / 3.0; 3];
        let result = soup_only(&toc, &profile, 0.05, 0.01);
        assert_eq!(result.rows(), 3);
        assert_eq!(result.cols(), 2);
        assert_eq!(result.nnz(), 0);
    }

    #[test]
    fn test_soup_only_zero_rho() {
        // With rho=0, no background → all entries kept
        let indptr = vec![0usize, 2, 3];
        let indices = vec![0usize, 1, 0];
        let data = vec![5.0, 3.0, 2.0];
        let toc = sprs::CsMatI::new((2, 2), indptr, indices, data);
        let profile = vec![0.5, 0.5];

        let result = soup_only(&toc, &profile, 0.0, 0.01);
        assert_eq!(result.nnz(), 3);
    }
}
