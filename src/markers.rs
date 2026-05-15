//! Quick marker gene discovery via tf-idf + hypergeometric test.
//!
//! Finds cluster-specific genes by:
//! 1. Binarizing the count matrix (counts > express_cut → expressed)
//! 2. Computing tf-idf scores per cluster × gene
//! 3. Hypergeometric test + BH correction, returning top N per cluster
//!
//! Parallelized across clusters with rayon.
//!
//! See REFACTOR_PLAN.md §1.4b for algorithm details.

use crate::stats::bh_correction;
use rayon::prelude::*;

/// Result for a single marker gene in a cluster.
#[derive(Debug, Clone)]
pub struct MarkerGene {
    pub gene_idx: usize,
    pub tfidf_score: f64,
    pub qval: f64,
}

/// TF-IDF based quick marker discovery.
///
/// # Arguments
/// * `toc` - Cell count matrix (CSR, genes × cells)
/// * `clusters` - Cluster assignment per cell (cluster IDs as usize)
/// * `n` - Number of top markers per cluster (default 10)
/// * `fdr` - False discovery rate threshold (default 0.01)
/// * `express_cut` - Binarization cutoff: counts > cut → expressed (default 0.9)
///
/// # Returns
/// Vec of marker lists, one per cluster. Each list is sorted by tfidf_score descending.
pub fn quick_markers(
    toc: &sprs::CsMatI<f64, usize>,
    clusters: &[usize],
    n: usize,
    fdr: f64,
    express_cut: f64,
) -> Vec<Vec<MarkerGene>> {
    let n_genes = toc.rows();
    let n_cells = clusters.len();
    let n_clusters = clusters.iter().max().map(|&m| m + 1).unwrap_or(0);

    if n_genes == 0 || n_cells == 0 || n_clusters == 0 {
        return vec![Vec::new(); n_clusters];
    }

    // ---------------------------------------------------------------
    // 1. Binarize: count > express_cut → gene is "expressed" in cell
    //    Build per-cluster expression counts
    // ---------------------------------------------------------------
    // gene_expr_total[gene] = total number of cells expressing this gene
    let mut gene_expr_total = vec![0usize; n_genes];

    // cluster_gene_expr[cluster][gene] = cells in cluster expressing gene
    let mut cluster_gene_expr: Vec<Vec<usize>> = vec![vec![0usize; n_genes]; n_clusters];
    let mut cluster_sizes = vec![0usize; n_clusters];

    // Iterate over non-zero entries
    for (&val, (gene_idx, col_idx)) in toc.iter() {
        if val > express_cut {
            let cluster = clusters.get(col_idx).copied().unwrap_or(0);
            gene_expr_total[gene_idx] += 1;
            if cluster < n_clusters {
                cluster_gene_expr[cluster][gene_idx] += 1;
            }
        }
    }

    // Count cluster sizes (cells per cluster, expressed or not)
    for &cluster in clusters {
        if cluster < n_clusters {
            cluster_sizes[cluster] += 1;
        }
    }

    // ---------------------------------------------------------------
    // 2. Compute tf-idf scores per cluster × gene (parallelized)
    // ---------------------------------------------------------------
    let results: Vec<Vec<MarkerGene>> = (0..n_clusters)
        .into_par_iter()
        .map(|cluster| {
            let cs = cluster_sizes[cluster] as f64;
            if cs == 0.0 {
                return Vec::new();
            }

            let genes: Vec<(usize, f64)> = (0..n_genes)
                .filter_map(|gene| {
                    let expressed = cluster_gene_expr[cluster][gene];
                    if expressed == 0 {
                        return None;
                    }
                    let tf = expressed as f64 / cs;
                    let total_expr = gene_expr_total[gene] as f64;
                    if total_expr <= 0.0 || total_expr >= n_cells as f64 {
                        return None;
                    }
                    let idf = (1.0 + (n_cells as f64) / total_expr).ln();
                    let score = tf * idf;
                    Some((gene, score))
                })
                .collect();

            // -----------------------------------------------------------
            // 3. Hypergeometric test + BH correction
            // -----------------------------------------------------------
            // H0: gene is not enriched in this cluster
            // P(X ≥ expressed) where X ~ Hypergeometric(N, K, n)
            //   N = total cells
            //   K = total cells expressing the gene
            //   n = cells in cluster
            //   k = cells in cluster expressing the gene
            let pvals: Vec<f64> = genes
                .iter()
                .map(|&(gene, _)| {
                    let k = cluster_gene_expr[cluster][gene];
                    let k_big = gene_expr_total[gene]; // K
                    let n = cluster_sizes[cluster]; // n (cluster size)
                    let big_n = n_cells; // N
                    hypergeometric_sf(k, k_big, n, big_n)
                })
                .collect();

            let qvals = bh_correction(&pvals);

            // Filter by FDR, keep top n by tfidf_score
            let mut markers: Vec<MarkerGene> = genes
                .into_iter()
                .enumerate()
                .filter(|(i, _)| qvals[*i] <= fdr)
                .map(|(i, (gene, score))| MarkerGene {
                    gene_idx: gene,
                    tfidf_score: score,
                    qval: qvals[i],
                })
                .collect();

            markers.sort_by(|a, b| b.tfidf_score.partial_cmp(&a.tfidf_score).unwrap_or(std::cmp::Ordering::Equal));
            markers.truncate(n);
            markers
        })
        .collect();

    results
}

/// Hypergeometric survival function: P(X ≥ k) under Hypergeometric(N, K, n).
///
/// Uses a numerically stable iterative computation avoiding large factorials.
/// For k ≤ K/2: uses CDF; for k > K/2: uses SF = 1 - CDF(k-1).
fn hypergeometric_sf(k: usize, k_big: usize, n: usize, big_n: usize) -> f64 {
    if k > k_big.min(n) || k == 0 {
        return 1.0;
    }

    // p(k) = C(K, k) * C(N-K, n-k) / C(N, n)
    // Compute iteratively from k to min(K, n)
    let max_k = k_big.min(n);

    // Compute p(k) using ln-space for stability, then convert
    let log_first = ln_choose(k_big, k) + ln_choose(big_n - k_big, n - k) - ln_choose(big_n, n);
    let p_k = log_first.exp();
    let mut sf = p_k;

    // Iterate k+1, k+2, ..., max_k
    let mut current = p_k;
    for i in k..max_k {
        // ratio = p(i+1) / p(i) = (K-i)*(n-i) / ((i+1)*(N-K-n+i+1))
        let num = (k_big - i) as f64 * (n - i) as f64;
        let den = (i + 1) as f64 * (big_n - k_big - n + i + 1) as f64;
        if den <= 0.0 {
            break;
        }
        current *= num / den;
        sf += current;
    }

    sf.min(1.0).max(0.0)
}

/// ln-choose: ln(C(n, k)) using the log-gamma function.
fn ln_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    // C(n, k) = n! / (k! * (n-k)!)
    // ln(C) = lnΓ(n+1) - lnΓ(k+1) - lnΓ(n-k+1)
    use statrs::function::gamma::ln_gamma;
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hypergeometric_sf_basic() {
        // P(X ≥ 5) in Hypergeometric(N=100, K=20, n=10)
        // Expected: small probability (5 is high draw from 20 out of 100)
        let p = hypergeometric_sf(5, 20, 10, 100);
        assert!(p > 0.0 && p < 0.1, "SF should be small: got {p}");
    }

    #[test]
    fn test_hypergeometric_sf_zero() {
        // P(X ≥ 0) = 1
        assert!((hypergeometric_sf(0, 10, 5, 100) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hypergeometric_sf_all() {
        // If K=n=big_N, then P(X ≥ n) = 1
        let p = hypergeometric_sf(5, 5, 5, 5);
        assert!((p - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ln_choose() {
        let val = ln_choose(10, 3);
        assert!(val > 0.0); // ln(C(10,3)) = ln(120) ≈ 4.79
        assert!((val - 4.787491742782046).abs() < 1e-6);
    }

    #[test]
    fn test_quick_markers_empty() {
        let toc = sprs::CsMatI::<f64, usize>::zero((5, 0));
        let clusters: Vec<usize> = vec![];
        let result = quick_markers(&toc, &clusters, 10, 0.01, 0.9);
        assert!(result.is_empty());
    }

    #[test]
    fn test_quick_markers_small() {
        // 3 genes × 6 cells, two clusters
        let indptr = vec![0usize, 3, 5, 8];
        let indices = vec![0usize, 2, 3, 0, 3, 1, 2, 4];
        let data = vec![1.0, 1.5, 4.0, 1.0, 2.0, 2.0, 1.0, 5.0];
        let toc = sprs::CsMatI::new((3, 6), indptr, indices, data);

        // 6 cells: 3 in cluster 0, 3 in cluster 1
        let clusters = vec![0usize, 0, 0, 1, 1, 1];

        let result = quick_markers(&toc, &clusters, 5, 1.0, 0.9);
        assert_eq!(result.len(), 2);
    }
}
