//! Contamination fraction (ρ) estimation.
//!
//! Two modes:
//! 1. Manual (`calculate_contamination_fraction`): user provides gene sets
//!    known to be absent in specific cell types → Poisson GLM → ρ
//! 2. Automatic (`auto_est_cont`): uses quickMarkers to find gene sets,
//!    then aggregates Gamma posterior densities across multiple estimates
//!    to find the consensus ρ (MAP in contamination_range).
//!
//! See REFACTOR_PLAN.md §1.4d-e for algorithm details.

use crate::markers::{self, MarkerGene};
use crate::non_expressing;
use crate::stats::gamma;

/// Manual contamination fraction estimation via Poisson GLM.
///
/// User provides gene sets known absent in some cell types.
/// For each gene set, identifies non-expressing cells, then
/// fits a Poisson regression to estimate ρ.
///
/// Returns: estimated ρ, or None if estimation fails.
pub fn calculate_contamination_fraction(
    toc: &sprs::CsMatI<f64, usize>,
    gene_sets: &[Vec<usize>],
    _clusters: &[usize],
    soup_profile: &[f64],
    max_contam: f64,
    fdr: f64,
) -> Option<f64> {
    if gene_sets.is_empty() {
        return Some(0.0);
    }

    let mut all_obs: Vec<f64> = Vec::new();
    let mut all_exp: Vec<f64> = Vec::new();

    for gene_indices in gene_sets {
        let non_expr = non_expressing::estimate_non_expressing_cells(
            toc, gene_indices, soup_profile, max_contam, fdr,
        );

        // Use non-expressing cells to estimate ρ
        // For these cells: obs ~ Poisson(exp * ρ) where exp = nUMIs * sum(soup[genes])
        let gene_set_soup: f64 = gene_indices
            .iter()
            .map(|&g| soup_profile.get(g).copied().unwrap_or(0.0))
            .sum();

        for (&val, (gene_idx, cell_idx)) in toc.iter() {
            if !non_expr.get(cell_idx).copied().unwrap_or(false) {
                continue;
            }
            if gene_indices.contains(&gene_idx) {
                let total = cell_total_umis(toc, cell_idx);
                all_obs.push(val);
                all_exp.push(total * gene_set_soup);
            }
        }
    }

    if all_obs.is_empty() || all_exp.is_empty() {
        return None;
    }

    // Simple Poisson rate estimation: ρ̂ = mean(obs) / mean(exp)
    let mean_obs: f64 = all_obs.iter().sum::<f64>() / all_obs.len() as f64;
    let mean_exp: f64 = all_exp.iter().sum::<f64>() / all_exp.len() as f64;

    if mean_exp <= 0.0 {
        return None;
    }

    let rho = (mean_obs / mean_exp).min(1.0).max(0.0);
    Some(rho)
}

/// Get total UMIs for a cell (sum of all gene counts).
fn cell_total_umis(toc: &sprs::CsMatI<f64, usize>, cell_idx: usize) -> f64 {
    let mut total = 0.0;
    for row in 0..toc.rows() {
        if let Some(val) = toc.get(row, cell_idx) {
            total += *val;
        }
    }
    total
}

/// Automatic contamination estimation with posterior density aggregation.
///
/// Algorithm:
/// 1. quickMarkers → marker genes per cluster
/// 2. For each marker gene × each cluster, estimate non-expressing cells
/// 3. Build Gamma posterior: Gamma(prior_k + obs, scale = prior_theta / (1 + prior_theta * exp))
/// 4. Aggregate all posteriors → find MAP in contamination_range
pub fn auto_est_cont(
    toc: &sprs::CsMatI<f64, usize>,
    clusters: &[usize],
    soup_profile: &[f64],
    tfidf_min: f64,
    _soup_quantile: f64,
    max_markers: usize,
    contamination_range: (f64, f64),
    prior_rho: f64,
    prior_rho_stddev: f64,
) -> Option<f64> {
    let (range_min, range_max) = contamination_range;
    let n_rho_probes = 1000;

    // -----------------------------------------------------------
    // 1. Find marker genes per cluster and build per-cluster gene GROUPS
    //    (using gene groups instead of single genes avoids Poisson test
    //    being too strict when individual soup proportions are ~1e-5)
    // -----------------------------------------------------------
    let all_markers: Vec<Vec<MarkerGene>> =
        markers::quick_markers(toc, clusters, max_markers, 0.05, 0.9);

    // Build per-cluster gene sets as Vec<Vec<usize>>
    let mut cluster_gene_sets: Vec<Vec<usize>> = Vec::new();
    for (_cluster, marker_list) in all_markers.iter().enumerate() {
        let mut gene_set: Vec<usize> = Vec::new();
        for m in marker_list {
            if m.tfidf_score >= tfidf_min {
                let soup_val = soup_profile.get(m.gene_idx).copied().unwrap_or(0.0);
                if soup_val > 0.0 {
                    gene_set.push(m.gene_idx);
                }
            }
        }
        cluster_gene_sets.push(gene_set);
    }

    // Adaptive fallback: lower tfidf_min if no clusters have enough genes
    let min_genes_per_cluster = 3usize;
    let total_genes: usize = cluster_gene_sets.iter().map(|gs| gs.len()).sum();
    if total_genes < 5 && tfidf_min > 0.1 {
        for fallback_min in [0.5f64, 0.2, 0.1] {
            if fallback_min >= tfidf_min {
                continue;
            }
            cluster_gene_sets.clear();
            for (_cluster, marker_list) in all_markers.iter().enumerate() {
                let mut gene_set: Vec<usize> = Vec::new();
                for m in marker_list {
                    if m.tfidf_score >= fallback_min {
                        let soup_val = soup_profile.get(m.gene_idx).copied().unwrap_or(0.0);
                        if soup_val > 0.0 {
                            gene_set.push(m.gene_idx);
                        }
                    }
                }
                cluster_gene_sets.push(gene_set);
            }
            let new_total: usize = cluster_gene_sets.iter().map(|gs| gs.len()).sum();
            if new_total >= 5 {
                break;
            }
        }
    }

    // Filter to non-empty gene sets
    let active_sets: Vec<&[usize]> = cluster_gene_sets
        .iter()
        .filter(|gs| gs.len() >= min_genes_per_cluster)
        .map(|gs| gs.as_slice())
        .collect();

    if active_sets.is_empty() {
        return Some(prior_rho);
    }

    // -----------------------------------------------------------
    // 2. For each active gene set, estimate non-expressing cells and accumulate
    // -----------------------------------------------------------
    let prior_theta = if prior_rho > 0.0 && prior_rho_stddev > 0.0 {
        prior_rho_stddev * prior_rho_stddev / prior_rho
    } else {
        1.0 / prior_rho.max(1e-6)
    };
    let prior_k = if prior_theta > 0.0 {
        prior_rho / prior_theta + 1.0
    } else {
        2.0
    };

    let mut posteriors: Vec<(f64, f64)> = Vec::new();

    // Precompute cell totals once
    let n_cells = toc.cols();
    let mut cell_totals = vec![0.0f64; n_cells];
    for (&val, (_g, cell_idx)) in toc.iter() {
        cell_totals[cell_idx] += val;
    }

    // Process each gene set as a group
    for gene_set in &active_sets {
        // Sum of soup proportions for this gene set
        let gene_set_soup: f64 = gene_set
            .iter()
            .map(|&g| soup_profile.get(g).copied().unwrap_or(0.0))
            .sum();

        if gene_set_soup <= 1e-12 {
            continue;
        }

        // Build O(1) membership mask
        let n_genes = toc.rows();
        let mut gene_mask = vec![false; n_genes];
        for &g in *gene_set {
            if g < n_genes { gene_mask[g] = true; }
        }

        // Find non-expressing cells for this gene set
        let non_expr = non_expressing::estimate_non_expressing_cells(
            toc, gene_set, soup_profile, range_max, 0.05,
        );

        // Accumulate obs and expected for non-expressing cells (single scan)
        let mut obs_sum = 0.0f64;
        let mut exp_sum = 0.0f64;

        for (&val, (g_idx, cell_idx)) in toc.iter() {
            if !non_expr.get(cell_idx).copied().unwrap_or(false) {
                continue;
            }
            if g_idx >= n_genes || !gene_mask[g_idx] {
                continue;
            }
            let total = cell_totals[cell_idx];
            obs_sum += val;
            exp_sum += total * gene_set_soup;
        }

        if exp_sum <= 0.0 || obs_sum <= 0.0 {
            continue;
        }

        // Gamma posterior: Gamma(prior_k + obs_sum, prior_theta / (1 + prior_theta * exp_sum))
        let post_k = prior_k + obs_sum;
        let post_theta = prior_theta / (1.0 + prior_theta * exp_sum);

        if post_theta > 0.0 {
            posteriors.push((post_k, post_theta));
        }
    }

    if posteriors.is_empty() {
        return Some(prior_rho);
    }

    // -----------------------------------------------------------
    // 3. Aggregate all Gamma posteriors into a combined density
    // -----------------------------------------------------------
    let mut best_rho = prior_rho;
    let mut best_log_density = f64::NEG_INFINITY;

    for i in 0..n_rho_probes {
        let rho = range_min + (range_max - range_min) * (i as f64) / ((n_rho_probes - 1) as f64);

        // Sum of log-densities from all posteriors (weighted equally)
        let total_log_density: f64 = posteriors
            .iter()
            .map(|&(k, theta)| gamma::log_dgamma(rho, k, theta))
            .sum();

        if total_log_density > best_log_density {
            best_log_density = total_log_density;
            best_rho = rho;
        }
    }

    // Clamp to valid range
    Some(best_rho.clamp(range_min, range_max))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_auto_est_cont_empty() {
        let toc = sprs::CsMatI::<f64, usize>::zero((3, 0));
        let clusters: Vec<usize> = vec![];
        let profile = vec![1.0 / 3.0; 3];
        let rho = auto_est_cont(
            &toc, &clusters, &profile, 1.0, 0.9, 50, (0.01, 0.8), 0.05, 0.10,
        );
        // With no cells, falls back to prior
        assert!(rho.is_some());
        assert!((rho.unwrap() - 0.05).abs() < 0.01);
    }

    #[test]
    fn test_calculate_contamination_empty_sets() {
        let toc = sprs::CsMatI::<f64, usize>::zero((5, 10));
        let clusters = vec![0usize; 10];
        let profile = vec![0.2; 5];
        let rho = calculate_contamination_fraction(
            &toc, &[], &clusters, &profile, 1.0, 0.05,
        );
        assert!(rho.is_some());
        assert!((rho.unwrap() - 0.0).abs() < 1e-10);
    }
}
