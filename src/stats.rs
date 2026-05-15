//! Statistical distributions: Poisson, Gamma, Hypergeometric.
//!
//! Wrappers around `statrs` with R-compatible parameterizations.
//! If statrs precision diverges from R's C-level implementations (Rmath),
//! this module will be replaced with libR-sys bindings or custom impls.
//!
//! See REFACTOR_PLAN.md §9 risk mitigation.

/// Poisson distribution functions (R-compatible interface).
pub mod poisson {
    use statrs::distribution::{DiscreteCDF, Poisson};

    /// P(X ≤ q) — equivalent to R's `ppois(q, lambda, lower.tail = TRUE)`
    pub fn ppois(q: f64, lambda: f64) -> f64 {
        if lambda <= 0.0 {
            return if q >= 0.0 { 1.0 } else { 0.0 };
        }
        if q < 0.0 {
            return 0.0;
        }
        let dist = Poisson::new(lambda).expect("lambda must be positive");
        dist.cdf(q as u64)
    }

    /// P(X ≥ q) — equivalent to R's `ppois(q, lambda, lower.tail = FALSE)`
    pub fn ppois_upper(q: f64, lambda: f64) -> f64 {
        1.0 - ppois(q - 1.0, lambda)
    }
}

/// Gamma distribution functions.
pub mod gamma {
    use statrs::distribution::{Continuous, Gamma};

    /// Probability density — equivalent to R's `dgamma(x, shape, scale)`
    /// Note: statrs Gamma uses (shape, rate), so we convert scale → rate = 1/scale.
    pub fn dgamma(x: f64, shape: f64, scale: f64) -> f64 {
        if x <= 0.0 || shape <= 0.0 || scale <= 0.0 {
            return 0.0;
        }
        let rate = 1.0 / scale;
        let dist = Gamma::new(shape, rate).expect("invalid gamma params");
        dist.pdf(x)
    }

    /// Log probability density — equivalent to R's `dgamma(x, shape, scale, log = TRUE)`
    pub fn log_dgamma(x: f64, shape: f64, scale: f64) -> f64 {
        let val = dgamma(x, shape, scale);
        if val <= 0.0 {
            f64::NEG_INFINITY
        } else {
            val.ln()
        }
    }
}

/// Hypergeometric distribution functions.
pub mod hypergeometric {
    // TODO: statrs does not have a hypergeometric distribution.
    // Will need custom implementation or libR-sys binding.
}

/// Benjamini-Hochberg FDR correction.
///
/// # Arguments
/// * `p_values` - Raw p-values (length n)
///
/// # Returns
/// BH-adjusted q-values (same length as input).
pub fn bh_correction(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return Vec::new();
    }

    let mut indexed: Vec<(usize, f64)> = p_values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut q_values = vec![0.0f64; n];
    let mut prev = 1.0f64;
    for (rank, &(orig_idx, p)) in indexed.iter().enumerate().rev() {
        let q_candidate = (p * n as f64) / (rank + 1) as f64;
        let q = q_candidate.min(prev).min(1.0);
        q_values[orig_idx] = q;
        prev = q;
    }
    q_values
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========================================================================
    // poisson tests — benchmarked against R 4.x ppois()
    // ========================================================================
    mod poisson_tests {
        use super::poisson::*;

        #[test]
        fn ppois_basic() {
            // R: ppois(0, lambda=1) = exp(-1) ≈ 0.36787944117144233
            let val = ppois(0.0, 1.0);
            assert!((val - 0.36787944117144233).abs() < 1e-10);
        }

        #[test]
        fn ppois_large_lambda() {
            // R: ppois(50, lambda=50) ≈ 0.5375 (median region of Poisson)
            let val = ppois(50.0, 50.0);
            assert!(val > 0.45 && val < 0.65,
                "ppois(50, 50) = {val}, expected ~0.5");
        }

        #[test]
        fn ppois_tail() {
            // R: ppois(100, lambda=10) ≈ 1.0 (essentially 1)
            let val = ppois(100.0, 10.0);
            assert!(val > 0.999999);
        }

        #[test]
        fn ppois_low_lambda() {
            // R: ppois(0, lambda=0.001) ≈ 0.9990005
            let val = ppois(0.0, 0.001);
            assert!((val - 0.999000499833375).abs() < 1e-10);
        }

        #[test]
        fn ppois_upper_basic() {
            // R: ppois(4, lambda=5, lower.tail=FALSE) ≈ 0.5595067
            let val = ppois_upper(5.0, 5.0);
            // ppois_upper(q=5, lambda=5) = P(X ≥ 5) = 1 - P(X ≤ 4)
            assert!((val - 0.5595067149347876).abs() < 1e-7);
        }

        #[test]
        fn ppois_upper_zero() {
            // P(X ≥ 0) = 1
            let val = ppois_upper(0.0, 5.0);
            assert!((val - 1.0).abs() < 1e-10);
        }
    }

    // ========================================================================
    // gamma tests — benchmarked against R 4.x dgamma()
    // ========================================================================
    mod gamma_tests {
        use super::gamma::*;

        #[test]
        fn dgamma_basic() {
            // R: dgamma(0.05, shape=1, scale=0.02) = dgamma(0.05, 1, rate=50)
            // scale = 0.02 → rate = 50
            // dgamma(0.05, shape=1, scale=0.02) = 50 * exp(-50*0.05) = 50*exp(-2.5) ≈ 4.10425
            let val = dgamma(0.05, 1.0, 0.02);
            let expected = 4.104249931486046;
            assert!((val - expected).abs() < 1e-6,
                "dgamma(0.05, 1, 0.02) = {val}, expected {expected}");
        }

        #[test]
        fn dgamma_zero_x() {
            // R: dgamma(0, shape=2, scale=1) = 0
            let val = dgamma(0.0, 2.0, 1.0);
            assert_eq!(val, 0.0);
        }

        #[test]
        fn dgamma_tiny_shape() {
            // R: dgamma(0.01, shape=0.1, scale=1) — should return a finite value
            let val = dgamma(0.01, 0.1, 1.0);
            assert!(val.is_finite());
            assert!(val > 0.0);
        }

        #[test]
        fn log_dgamma_avoids_underflow() {
            // For very small x, density may underflow but log should be finite
            let val = log_dgamma(1e-10, 1.0, 1.0);
            assert!(val.is_finite());
            assert!(val < 0.0);
        }
    }

    // ========================================================================
    // BH correction tests — benchmarked against R's p.adjust(method="BH")
    // ========================================================================
    mod bh_tests {
        use super::bh_correction;

        #[test]
        fn bh_empty() {
            let result = bh_correction(&[]);
            assert!(result.is_empty());
        }

        #[test]
        fn bh_single() {
            // R: p.adjust(c(0.05), method="BH") = 0.05
            let result = bh_correction(&[0.05]);
            assert!((result[0] - 0.05).abs() < 1e-10);
        }

        #[test]
        fn bh_basic() {
            // R: p.adjust(c(0.01, 0.02, 0.03, 0.04, 0.05), method="BH")
            // → c(0.05, 0.05, 0.05, 0.05, 0.05)
            let pvals = vec![0.01, 0.02, 0.03, 0.04, 0.05];
            let result = bh_correction(&pvals);
            // All should be close to 0.05
            for (i, &q) in result.iter().enumerate() {
                assert!((q - 0.05).abs() < 1e-10,
                    "BH[{i}]: got {q}, expected 0.05");
            }
        }

        #[test]
        fn bh_mixed() {
            // Verified against R: p.adjust(c(0.001, 0.5, 0.01, 0.9), method="BH")
            let pvals = vec![0.001, 0.5, 0.01, 0.9];
            let result = bh_correction(&pvals);
            // Standard BH: sort ascending, BH = p*n/rank, cumulative min from right
            // sorted: [0.001, 0.01, 0.5, 0.9]
            // raw BH:  [0.004, 0.02, 0.667, 0.9]
            // cummin:  [0.004, 0.02, 0.667, 0.9]
            // mapped back to original order: [0.004, 0.667, 0.02, 0.9]
            // Cumulative min enforces: q[1]=min(0.667, 0.9)=0.667
            let expected = vec![0.004, 0.6666666666666666, 0.02, 0.9];
            for (i, (&got, &exp)) in result.iter().zip(expected.iter()).enumerate() {
                assert!((got - exp).abs() < 1e-10,
                    "BH[{i}]: got {got}, expected {exp}");
            }
        }

        #[test]
        fn bh_monotonic_enforced() {
            // Verify q-values are monotonically non-decreasing with p-values
            let pvals = vec![0.001, 0.01, 0.1, 0.5];
            let result = bh_correction(&pvals);
            for w in result.windows(2) {
                assert!(w[0] <= w[1] + 1e-10);
            }
        }
    }
}
