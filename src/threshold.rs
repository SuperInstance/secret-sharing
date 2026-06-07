//! Threshold secret sharing utilities.
//!
//! Higher-level interface for threshold-based secret sharing with
//! convenience methods for common configurations.

use crate::field::FiniteField;
use crate::shamir::{ShamirScheme, Share};
use crate::reconstruct::Reconstructor;

/// Threshold secret sharing configuration.
#[derive(Debug, Clone)]
pub struct ThresholdScheme {
    /// Number of total shares.
    pub n: usize,
    /// Threshold for reconstruction.
    pub threshold: usize,
}

impl ThresholdScheme {
    /// Create a new threshold scheme with n shares and threshold t.
    pub fn new(n: usize, t: usize) -> Self {
        assert!(n >= t, "n must be >= t");
        assert!(t >= 2, "threshold must be >= 2 for security");
        Self { n, threshold: t }
    }

    /// Create a (t, n) = (n/2, n) scheme (majority threshold).
    pub fn majority(n: usize) -> Self {
        Self::new(n, (n + 1).div_ceil(2))
    }

    /// Create a (t, n) = (2, n) scheme (any 2 of n).
    pub fn any_two(n: usize) -> Self {
        Self::new(n, 2)
    }

    /// Split a secret using this threshold scheme.
    pub fn split(&self, secret: &FiniteField) -> Vec<Share> {
        ShamirScheme::split(secret, self.n, self.threshold)
    }

    /// Reconstruct from a set of shares.
    pub fn reconstruct(&self, shares: &[Share]) -> FiniteField {
        assert!(shares.len() >= self.threshold, "not enough shares");
        Reconstructor::reconstruct(shares)
    }

    /// Verify that the given number of shares is sufficient.
    pub fn is_sufficient(&self, share_count: usize) -> bool {
        share_count >= self.threshold
    }

    /// Get the security level (threshold / total shares ratio).
    pub fn security_ratio(&self) -> f64 {
        self.threshold as f64 / self.n as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn majority_scheme() {
        let scheme = ThresholdScheme::majority(5);
        assert_eq!(scheme.threshold, 3);
        assert_eq!(scheme.n, 5);
    }

    #[test]
    fn any_two_scheme() {
        let scheme = ThresholdScheme::any_two(6);
        assert_eq!(scheme.threshold, 2);
        assert_eq!(scheme.n, 6);
    }

    #[test]
    fn split_and_reconstruct_majority() {
        let scheme = ThresholdScheme::majority(5);
        let secret = FiniteField::from_u64(42);
        let shares = scheme.split(&secret);
        assert_eq!(shares.len(), 5);

        let subset: Vec<Share> = vec![shares[0].clone(), shares[2].clone(), shares[4].clone()];
        let recovered = scheme.reconstruct(&subset);
        assert_eq!(recovered, secret);
    }

    #[test]
    fn is_sufficient() {
        let scheme = ThresholdScheme::new(5, 3);
        assert!(scheme.is_sufficient(3));
        assert!(scheme.is_sufficient(5));
        assert!(!scheme.is_sufficient(2));
    }

    #[test]
    #[should_panic(expected = "threshold must be >= 2")]
    fn threshold_too_low() {
        ThresholdScheme::new(5, 1);
    }

    #[test]
    #[should_panic(expected = "not enough shares")]
    fn reconstruct_insufficient_shares() {
        let scheme = ThresholdScheme::new(5, 3);
        let secret = FiniteField::from_u64(42);
        let shares = scheme.split(&secret);
        let insufficient: Vec<Share> = shares.into_iter().take(2).collect();
        scheme.reconstruct(&insufficient);
    }
}
