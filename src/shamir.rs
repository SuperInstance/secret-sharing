//! Shamir's secret sharing scheme.
//!
//! Split a secret into `n` shares such that any `threshold` shares can
//! reconstruct it, but fewer than `threshold` reveal nothing.

use crate::field::FiniteField;
use num_bigint::BigInt;
use rand::Rng;

/// A single share in Shamir's scheme: (x, y) where x is the share index.
#[derive(Debug, Clone, PartialEq)]
pub struct Share {
    /// The x-coordinate (share index, typically 1..=n).
    pub x: BigInt,
    /// The y-coordinate (the share value).
    pub y: FiniteField,
}

/// Shamir's secret sharing scheme.
pub struct ShamirScheme;

impl ShamirScheme {
    /// Split a secret into `n` shares with threshold `t`.
    ///
    /// Any `t` shares can reconstruct the secret. `n` must be >= `t`.
    pub fn split(secret: &FiniteField, n: usize, t: usize) -> Vec<Share> {
        assert!(n >= t, "n must be >= threshold t");
        assert!(t >= 1, "threshold must be >= 1");

        let modulus = secret.modulus().clone();

        // Generate random polynomial coefficients: f(x) = secret + a1*x + a2*x^2 + ... + a_{t-1}*x^{t-1}
        let mut coeffs = vec![secret.clone()];
        let mut rng = rand::thread_rng();
        for _ in 1..t {
            let rand_val: u64 = rng.gen();
            coeffs.push(FiniteField::new(BigInt::from(rand_val), modulus.clone()));
        }

        // Evaluate polynomial at x = 1, 2, ..., n
        let mut shares = Vec::with_capacity(n);
        for i in 1..=n {
            let x = BigInt::from(i);
            let y = Self::eval_poly(&coeffs, &x, &modulus);
            shares.push(Share { x: x.clone(), y });
        }
        shares
    }

    /// Evaluate a polynomial at point x in the finite field.
    fn eval_poly(coeffs: &[FiniteField], x: &BigInt, modulus: &BigInt) -> FiniteField {
        let mut result = FiniteField::zero(modulus);
        let mut x_power = FiniteField::one(modulus);
        for coeff in coeffs {
            let term = coeff.mul(&x_power);
            result = result.add(&term);
            x_power = x_power.mul(&FiniteField::new(x.clone(), modulus.clone()));
        }
        result
    }

    /// Compute the threshold for a given share count.
    pub fn min_threshold(n: usize) -> usize {
        (n + 1).div_ceil(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_produces_correct_count() {
        let secret = FiniteField::from_u64(12345);
        let shares = ShamirScheme::split(&secret, 5, 3);
        assert_eq!(shares.len(), 5);
    }

    #[test]
    fn split_threshold_1() {
        let secret = FiniteField::from_u64(42);
        let shares = ShamirScheme::split(&secret, 3, 1);
        assert_eq!(shares.len(), 3);
        // With threshold 1, every share is the secret itself
        for share in &shares {
            assert_eq!(share.y, secret);
        }
    }

    #[test]
    fn split_deterministic_secret_recovery() {
        let secret = FiniteField::from_u64(999);
        let shares = ShamirScheme::split(&secret, 5, 3);
        // Use first 3 shares to reconstruct
        let subset: Vec<&Share> = shares.iter().take(3).collect();
        let recovered = crate::reconstruct::Reconstructor::lagrange_interpolate(&subset);
        assert_eq!(recovered, secret);
    }

    #[test]
    fn any_subset_recovers() {
        let secret = FiniteField::from_u64(777);
        let shares = ShamirScheme::split(&secret, 5, 3);
        // Try different subsets
        let subsets: Vec<Vec<usize>> = vec![
            vec![0, 1, 2],
            vec![2, 3, 4],
            vec![0, 2, 4],
            vec![1, 3, 4],
        ];
        for indices in &subsets {
            let subset: Vec<&Share> = indices.iter().map(|&i| &shares[i]).collect();
            let recovered = crate::reconstruct::Reconstructor::lagrange_interpolate(&subset);
            assert_eq!(recovered, secret, "failed for subset {:?}", indices);
        }
    }

    #[test]
    #[should_panic(expected = "n must be >= threshold t")]
    fn split_panics_n_less_than_t() {
        let secret = FiniteField::from_u64(42);
        ShamirScheme::split(&secret, 2, 5);
    }

    #[test]
    fn large_threshold() {
        let secret = FiniteField::from_u64(1234);
        let shares = ShamirScheme::split(&secret, 10, 10);
        let subset: Vec<&Share> = shares.iter().collect();
        let recovered = crate::reconstruct::Reconstructor::lagrange_interpolate(&subset);
        assert_eq!(recovered, secret);
    }
}
