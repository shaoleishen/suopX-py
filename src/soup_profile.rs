//! Soup profile estimation.
//!
//! Estimates the proportion of each gene in ambient RNA ("soup")
//! from empty droplets (UMI range [0, 100]).
//!
//! This is deliberately kept simple and may remain a Python-level operation
//! since it only requires summation over a subset of droplets.

/// Estimate soup profile from empty droplet counts.
///
/// # Arguments
/// * `tod` - Full droplet matrix (CSR, genes × droplets)
/// * `empty_mask` - Boolean mask indicating empty droplets
///
/// # Returns
/// * `est` - Proportion of each gene in the soup (length = n_genes)
/// * `counts` - Raw UMI counts per gene in empty droplets
pub fn estimate_soup(
    _tod: &sprs::CsMatI<f64, usize>,
    _empty_mask: &[bool],
) -> (Vec<f64>, Vec<f64>) {
    // Placeholder - TODO: implement from R's estimateSoup
    unimplemented!("soup profile estimation")
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_empty_input() {
        // TODO
    }
}
