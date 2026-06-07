//! Lagrange interpolation for secret reconstruction.
//!
//! Given a sufficient number of shares, reconstruct the original secret
//! using Lagrange basis polynomials evaluated at x=0.

use crate::field::FiniteField;
use crate::shamir::Share;
use num_bigint::BigInt;

/// Secret reconstructor using Lagrange interpolation.
pub struct Reconstructor;

impl Reconstructor {
    /// Reconstruct the secret from shares using Lagrange interpolation at x=0.
    ///
    /// The number of shares must be at least the threshold used during splitting.
    pub fn lagrange_interpolate(shares: &[&Share]) -> FiniteField {
        assert!(!shares.is_empty(), "need at least one share");

        let modulus = shares[0].y.modulus().clone();
        let mut result = FiniteField::zero(&modulus);

        for (i, si) in shares.iter().enumerate() {
            let mut numerator = FiniteField::one(&modulus);
            let mut denominator = FiniteField::one(&modulus);

            for (j, sj) in shares.iter().enumerate() {
                if i == j {
                    continue;
                }
                // numerator *= (0 - x_j) = -x_j
                let neg_xj = FiniteField::new(BigInt::from(0) - &sj.x, modulus.clone());
                numerator = numerator.mul(&neg_xj);

                // denominator *= (x_i - x_j)
                let diff = FiniteField::new(&si.x - &sj.x, modulus.clone());
                denominator = denominator.mul(&diff);
            }

            let lagrange_coeff = numerator.div(&denominator).expect("denominator should not be zero");
            let term = si.y.mul(&lagrange_coeff);
            result = result.add(&term);
        }

        result
    }

    /// Reconstruct from owned shares.
    pub fn reconstruct(shares: &[Share]) -> FiniteField {
        let refs: Vec<&Share> = shares.iter().collect();
        Self::lagrange_interpolate(&refs)
    }

    /// Verify that a given set of shares reconstructs to an expected secret.
    pub fn verify(shares: &[&Share], expected: &FiniteField) -> bool {
        let reconstructed = Self::lagrange_interpolate(shares);
        reconstructed == *expected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shamir::ShamirScheme;

    #[test]
    fn reconstruct_from_minimum_shares() {
        let secret = FiniteField::from_u64(42);
        let shares = ShamirScheme::split(&secret, 5, 3);
        let subset: Vec<&Share> = shares.iter().take(3).collect();
        let recovered = Reconstructor::lagrange_interpolate(&subset);
        assert_eq!(recovered, secret);
    }

    #[test]
    fn reconstruct_from_all_shares() {
        let secret = FiniteField::from_u64(9999);
        let shares = ShamirScheme::split(&secret, 4, 2);
        let refs: Vec<&Share> = shares.iter().collect();
        let recovered = Reconstructor::lagrange_interpolate(&refs);
        assert_eq!(recovered, secret);
    }

    #[test]
    fn verify_correct() {
        let secret = FiniteField::from_u64(42);
        let shares = ShamirScheme::split(&secret, 5, 3);
        let subset: Vec<&Share> = shares.iter().take(3).collect();
        assert!(Reconstructor::verify(&subset, &secret));
    }

    #[test]
    fn verify_wrong_secret() {
        let secret = FiniteField::from_u64(42);
        let shares = ShamirScheme::split(&secret, 5, 3);
        let subset: Vec<&Share> = shares.iter().take(3).collect();
        let wrong = FiniteField::from_u64(99);
        assert!(!Reconstructor::verify(&subset, &wrong));
    }

    #[test]
    fn reconstruct_two_shares_threshold_two() {
        let secret = FiniteField::from_u64(1234);
        let shares = ShamirScheme::split(&secret, 4, 2);
        let subset: Vec<&Share> = vec![&shares[0], &shares[3]];
        let recovered = Reconstructor::lagrange_interpolate(&subset);
        assert_eq!(recovered, secret);
    }
}
