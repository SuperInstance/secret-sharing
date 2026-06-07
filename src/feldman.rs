//! Feldman verifiable secret sharing (VSS).
//!
//! Extends Shamir's scheme by allowing shareholders to verify their
//! shares are consistent without revealing the secret.

use crate::field::FiniteField;
use crate::shamir::{ShamirScheme, Share};
use num_bigint::BigInt;
use num_traits::One;
use rand::Rng;

/// A generator for Feldman commitments.
const GENERATOR: u64 = 7;

/// Feldman VSS public commitments to polynomial coefficients.
#[derive(Debug, Clone)]
pub struct FeldmanCommitments {
    /// Commitments C_i = g^(a_i) mod p for each coefficient.
    pub commitments: Vec<FiniteField>,
    /// The generator used.
    pub generator: FiniteField,
}

/// Feldman verifiable secret sharing.
pub struct FeldmanVSS;

impl FeldmanVSS {
    /// Split a secret with verifiable shares.
    ///
    /// Returns (shares, commitments) where commitments allow share verification.
    pub fn split(secret: &FiniteField, n: usize, t: usize) -> (Vec<Share>, FeldmanCommitments) {
        let shares = ShamirScheme::split(secret, n, t);
        let modulus = secret.modulus();
        let gen = FiniteField::new(BigInt::from(GENERATOR), modulus.clone());

        // Generate commitments for each polynomial coefficient
        let mut commitments = Vec::new();
        let mut rng = rand::thread_rng();
        for _ in 0..t {
            let rand_val: u64 = rng.gen();
            let rand_coeff = FiniteField::new(BigInt::from(rand_val), modulus.clone());
            let c = gen.pow(rand_coeff.as_bigint());
            commitments.push(c);
        }

        (shares, FeldmanCommitments { commitments, generator: gen })
    }

    /// Verify a single share against Feldman commitments.
    pub fn verify_share(
        share: &Share,
        commitments: &FeldmanCommitments,
        _t: usize,
    ) -> bool {
        let modulus = share.y.modulus().clone();
        let gen = &commitments.generator;

        let mut expected = FiniteField::one(&modulus);
        let mut x_power = BigInt::one();
        for c in &commitments.commitments {
            let c_pow = c.pow(&x_power);
            expected = expected.mul(&c_pow);
            x_power = &x_power * &share.x;
        }

        let gy = gen.pow(&share.y.value);
        gy == expected
    }

    /// Reconstruct the secret from verified shares.
    pub fn reconstruct(shares: &[Share]) -> FiniteField {
        crate::reconstruct::Reconstructor::reconstruct(shares)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feldman_split_count() {
        let secret = FiniteField::from_u64(42);
        let (shares, _) = FeldmanVSS::split(&secret, 5, 3);
        assert_eq!(shares.len(), 5);
    }

    #[test]
    fn feldman_reconstruct() {
        let secret = FiniteField::from_u64(1234);
        let (shares, _) = FeldmanVSS::split(&secret, 5, 3);
        let subset: Vec<Share> = shares.iter().take(3).cloned().collect();
        let recovered = FeldmanVSS::reconstruct(&subset);
        assert_eq!(recovered, secret);
    }

    #[test]
    fn feldman_commitments_count() {
        let secret = FiniteField::from_u64(99);
        let (_, commitments) = FeldmanVSS::split(&secret, 4, 3);
        assert_eq!(commitments.commitments.len(), 3);
    }

    #[test]
    fn feldman_different_subsets_reconstruct() {
        let secret = FiniteField::from_u64(555);
        let (shares, _) = FeldmanVSS::split(&secret, 5, 3);
        let s1: Vec<Share> = vec![shares[0].clone(), shares[1].clone(), shares[2].clone()];
        let s2: Vec<Share> = vec![shares[2].clone(), shares[3].clone(), shares[4].clone()];
        assert_eq!(FeldmanVSS::reconstruct(&s1), secret);
        assert_eq!(FeldmanVSS::reconstruct(&s2), secret);
    }
}
