//! Estimate non-expressing cells for a given gene set.
//!
//! Given a set of genes (e.g., HB genes expected absent in non-erythroid cells),
//! identifies cells that do not express them, using Poisson test with BH correction.
//!
//! Used as input to contamination fraction estimation.
//!
//! See REFACTOR_PLAN.md §1.4c for algorithm details.

use crate::stats::{bh_correction, poisson};

/// Identify cells that do not express a given gene set.
///
/// For each cell, computes expected background counts from soup profile,
/// then tests whether observed counts are consistent with pure contamination.
///
/// # Arguments
/// * `toc` - Cell count matrix (CSR, genes × cells)
/// * `gene_indices` - Indices of genes in the test set (global gene indices)
/// * `soup_profile` - Estimated soup proportions (length = n_genes)
/// * `max_contam` - Maximum plausible contamination fraction (upper bound for expected)
/// * `fdr` - False discovery rate threshold for BH correction
///
/// # Returns
/// Boolean mask: true if cell is determined to NOT express the gene set.
pub fn estimate_non_expressing_cells(
    toc: &sprs::CsMatI<f64, usize>,
    gene_indices: &[usize],
    soup_profile: &[f64],
    max_contam: f64,
    fdr: f64,
) -> Vec<bool> {
    let n_cells = toc.cols();

    if gene_indices.is_empty() || n_cells == 0 {
        return vec![false; n_cells];
    }

    // Sum of soup proportions for the gene set
    let gene_set_soup: f64 = gene_indices
        .iter()
        .map(|&g| soup_profile.get(g).copied().unwrap_or(0.0))
        .sum();

    if gene_set_soup <= 0.0 {
        return vec![true; n_cells]; // No background → all cells "non-expressing"
    }

    // For each cell: observed UMIs for gene set, total UMIs
    // Poisson test: P(X ≥ obs-1 | λ = total_UMIs * max_contam * gene_set_soup)
    let mut pvals = vec![1.0f64; n_cells];
    let mut cell_totals = vec![0.0f64; n_cells];
    let mut cell_gene_set = vec![0.0f64; n_cells];

    // Extract per-cell totals and gene set counts from CSR
    // Use Vec<bool> mask for O(1) gene set membership check
    let n_genes = toc.rows();
    let mut gene_set_mask = vec![false; n_genes];
    for &g in gene_indices {
        if g < n_genes {
            gene_set_mask[g] = true;
        }
    }

    for (&val, (gene_idx, cell_idx)) in toc.iter() {
        cell_totals[cell_idx] += val;
        if gene_idx < n_genes && gene_set_mask[gene_idx] {
            cell_gene_set[cell_idx] += val;
        }
    }

    for (cell_idx, (&total, &gs_counts)) in
        cell_totals.iter().zip(cell_gene_set.iter()).enumerate()
    {
        if total <= 0.0 {
            pvals[cell_idx] = 1.0;
            continue;
        }
        let expected = total * max_contam * gene_set_soup;
        if expected <= 0.0 {
            pvals[cell_idx] = 0.0; // Unexpected counts with zero expectation → significant
            continue;
        }
        // P(X ≥ max(obs-1, 0) | Poisson(expected))
        let obs_adj = (gs_counts - 1.0).max(0.0);
        pvals[cell_idx] = poisson::ppois_upper(obs_adj, expected);
    }

    // BH correction
    let qvals = bh_correction(&pvals);

    // Non-expressing if qval > fdr (fail to reject H0: counts are from background)
    qvals.iter().map(|&q| q > fdr).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_gene_set() {
        let toc = sprs::CsMatI::<f64, usize>::zero((3, 10));
        let result = estimate_non_expressing_cells(
            &toc, &[], &[0.1, 0.2, 0.3], 1.0, 0.05,
        );
        assert_eq!(result.len(), 10);
        assert!(result.iter().all(|&x| !x)); // Empty gene set → no decision
    }

    #[test]
    fn test_no_observed_counts() {
        // All cells have zero counts for the gene set → all non-expressing
        let toc = sprs::CsMatI::<f64, usize>::zero((5, 3));
        let profile = vec![1.0 / 5.0; 5];
        let result = estimate_non_expressing_cells(
            &toc, &[0, 1], &profile, 0.5, 0.05,
        );
        // With no counts, p=1.0 → q=1.0 → all non-expressing (q > 0.05)
        assert!(result.iter().all(|&x| x));
    }

    #[test]
    fn test_high_expression() {
        // 1 gene × 2 cells, both cells have high counts for this gene
        let indptr = vec![0usize, 2];
        let indices = vec![0usize, 1];
        let data = vec![100.0, 50.0];
        let toc = sprs::CsMatI::new((1, 2), indptr, indices, data);
        let profile = vec![0.01];
        let result = estimate_non_expressing_cells(
            &toc, &[0], &profile, 0.5, 0.05,
        );
        // Expected for cell 0: 100 * 0.5 * 0.01 = 0.5, obs=100 → super significant → EXPRESSING
        assert_eq!(result, vec![false, false]);
    }
}
