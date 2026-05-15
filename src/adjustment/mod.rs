//! Count adjustment engine — entry point.
//!
//! Three methods (see REFACTOR_PLAN.md §1.4f):
//! - subtraction: iterative per-gene background subtraction
//! - soup_only: Poisson test to identify and remove pure-contamination genes
//! - multinomial: explicit multinomial likelihood maximization (greedy)

pub mod subtraction;
pub mod soup_only;
pub mod multinomial;

/// Adjust counts to remove ambient RNA contamination.
///
/// # Arguments
/// * `toc` - Cell count matrix (CSR, genes × cells)
/// * `soup_profile` - Estimated soup proportions per gene
/// * `rho` - Global contamination fraction
/// * `method` - "subtraction", "soup_only", or "multinomial"
/// * `p_cut` - p-value cutoff (only used by soup_only)
///
/// # Returns
/// Corrected count matrix (CSR) with same sparsity pattern.
pub fn adjust_counts(
    toc: &sprs::CsMatI<f64, usize>,
    soup_profile: &[f64],
    rho: f64,
    method: &str,
    p_cut: f64,
) -> sprs::CsMatI<f64, usize> {
    match method {
        "subtraction" => subtraction::subtract(toc, soup_profile, rho),
        "soup_only" => soup_only::soup_only(toc, soup_profile, rho, p_cut),
        "multinomial" => multinomial::multinomial(toc, soup_profile, rho),
        _ => panic!("Unknown adjustment method: {method}"),
    }
}
