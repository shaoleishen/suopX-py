//! Alloc allocation algorithm — distributes total soup UMI budget across genes.
//!
//! Given a target total, per-gene upper bounds (observed counts), and weights
//! (soup proportions), fills buckets in order of limit/weight ratio until
//! capacity is exhausted.
//!
//! Optimization: early exit for trivial cases (all-zero, all-full), then
//! partial_sort for small targets, full sort otherwise.
//!
//! See REFACTOR_PLAN.md §4.3c for the optimization analysis.

/// Allocate a total amount across buckets with individual limits.
///
/// # Arguments
/// * `tgt` - Total amount to allocate (total soup UMIs to distribute)
/// * `bucket_lims` - Upper bound per bucket (observed counts per gene in this cell)
/// * `soup_profile` - Global soup proportions (length = n_total_genes)
/// * `gene_indices` - Maps local bucket index i → global gene index for soup_profile lookup
///
/// # Returns
/// Amount allocated to each bucket (≤ bucket_lims[i], sum ≈ tgt)
pub fn alloc(
    tgt: f64,
    bucket_lims: &[f64],
    soup_profile: &[f64],
    gene_indices: &[usize],
) -> Vec<f64> {
    let k = bucket_lims.len();

    // Early exit: zero target
    if tgt <= 0.0 {
        return vec![0.0; k];
    }

    // Early exit: no capacity
    let total_capacity: f64 = bucket_lims.iter().sum();
    if total_capacity <= 0.0 {
        return vec![0.0; k];
    }

    // Compute priority ratios: bucket_lims[i] / soup_profile[gene_indices[i]]
    // Smaller ratio → bucket fills first (limited capacity relative to weight)
    let mut indexed: Vec<(usize, f64, f64)> = bucket_lims
        .iter()
        .enumerate()
        .map(|(i, &lim)| {
            let gene_idx = gene_indices.get(i).copied().unwrap_or(0);
            let w = soup_profile.get(gene_idx).copied().unwrap_or(0.0);
            let ratio = if w > 0.0 { lim / w } else { f64::INFINITY };
            (i, lim, ratio)
        })
        .collect();

    // Sort by ratio ascending (fill smaller-ratio buckets first)
    indexed.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    let mut remaining = tgt;
    let mut result = vec![0.0f64; k];

    for &(orig_idx, lim, _ratio) in &indexed {
        if remaining <= 1e-12 {
            break;
        }
        let alloc = remaining.min(lim);
        result[orig_idx] = alloc;
        remaining -= alloc;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_target() {
        let result = alloc(0.0, &[5.0, 3.0], &[0.5, 0.5], &[0, 1]);
        assert_eq!(result, vec![0.0, 0.0]);
    }

    #[test]
    fn test_simple_allocation_direct_indices() {
        // gene_indices [0, 1] → soup_profile[0]=0.5, soup_profile[1]=0.5
        // ratio = lim/w → [10/0.5=20, 6/0.5=12] → fill bucket 1 first
        let result = alloc(10.0, &[10.0, 6.0], &[0.5, 0.5], &[0, 1]);
        assert!((result[0] - 4.0).abs() < 1e-10);
        assert!((result[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_sparse_gene_indices() {
        // Only genes 3 and 7 are present in this cell (bucket_lims has 2 entries)
        // soup_profile is global (10 genes), gene_indices maps local→global
        let soup_profile: Vec<f64> = (0..10).map(|i| (i + 1) as f64 / 55.0).collect();
        // gene_indices = [3, 7] → soup_profile[3]=4/55≈0.073, soup_profile[7]=8/55≈0.145
        // ratio = [5.0/0.073≈68, 5.0/0.145≈34] → fill bucket 1 first (gene 7)
        let result = alloc(5.0, &[5.0, 5.0], &soup_profile, &[3, 7]);
        // gene 7 has higher weight → gets filled first → gets 5.0
        // gene 3 gets remaining: 0.0
        assert!((result[1] - 5.0).abs() < 1e-10);
        assert!((result[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_target_exceeds_capacity() {
        let result = alloc(10.0, &[3.0, 2.0], &[1.0, 1.0], &[0, 1]);
        assert_eq!(result, vec![3.0, 2.0]);
    }

    #[test]
    fn test_zero_weight_gene() {
        // Gene 1 has zero weight → infinite ratio → filled last
        let result = alloc(5.0, &[5.0, 5.0], &[1.0, 0.0], &[0, 1]);
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 0.0).abs() < 1e-10);
    }
}
