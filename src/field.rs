//! Finite field arithmetic using modular arithmetic over a prime modulus.
//!
//! Provides a [`FiniteField`] type supporting addition, subtraction,
//! multiplication, inversion, and division modulo a prime.

use num_bigint::BigInt;
use num_traits::{One, Zero, ToPrimitive};
use std::fmt;

/// A large prime commonly used for secret sharing (Mersenne-like safe prime).
pub const DEFAULT_PRIME: &str = "170141183460469231731687303715884105727";

/// Finite field element with modular arithmetic.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FiniteField {
    /// The value (always in [0, modulus)).
    pub value: BigInt,
    /// The prime modulus.
    pub modulus: BigInt,
}

impl FiniteField {
    /// Create a new field element, reducing modulo the given prime.
    pub fn new(value: BigInt, modulus: BigInt) -> Self {
        let v = if value < BigInt::zero() {
            ((value % &modulus) + &modulus) % &modulus
        } else {
            value % &modulus
        };
        Self { value: v, modulus }
    }

    /// Create using the default prime.
    pub fn from_u64(value: u64) -> Self {
        let modulus = DEFAULT_PRIME.parse::<BigInt>().unwrap();
        Self::new(BigInt::from(value), modulus)
    }

    /// Create with a custom prime modulus.
    pub fn from_u64_with_modulus(value: u64, modulus: BigInt) -> Self {
        Self::new(BigInt::from(value), modulus)
    }

    /// Field element zero.
    pub fn zero(modulus: &BigInt) -> Self {
        Self { value: BigInt::zero(), modulus: modulus.clone() }
    }

    /// Field element one.
    pub fn one(modulus: &BigInt) -> Self {
        Self { value: BigInt::one(), modulus: modulus.clone() }
    }

    /// Get the modulus.
    pub fn modulus(&self) -> &BigInt {
        &self.modulus
    }

    /// Add two field elements.
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "modulus mismatch");
        Self::new(&self.value + &other.value, self.modulus.clone())
    }

    /// Subtract two field elements.
    pub fn sub(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "modulus mismatch");
        Self::new(&self.value - &other.value, self.modulus.clone())
    }

    /// Multiply two field elements.
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "modulus mismatch");
        Self::new(&self.value * &other.value, self.modulus.clone())
    }

    /// Compute the multiplicative inverse using extended GCD.
    pub fn inv(&self) -> Option<Self> {
        if self.value.is_zero() {
            return None;
        }
        // Extended Euclidean algorithm
        let (mut old_r, mut r) = (self.modulus.clone(), self.value.clone());
        let (mut old_s, mut s) = (BigInt::zero(), BigInt::one());

        while !r.is_zero() {
            let q = &old_r / &r;
            let new_r = &old_r - &q * &r;
            old_r = r;
            r = new_r;
            let new_s = &old_s - &q * &s;
            old_s = s;
            s = new_s;
        }

        if old_r != BigInt::one() {
            return None;
        }
        Some(Self::new(old_s, self.modulus.clone()))
    }

    /// Divide two field elements (multiply by inverse).
    pub fn div(&self, other: &Self) -> Option<Self> {
        other.inv().map(|inv| self.mul(&inv))
    }

    /// Exponentiate in the field.
    pub fn pow(&self, exp: &BigInt) -> Self {
        if exp.is_zero() {
            return Self::one(&self.modulus);
        }
        if exp < &BigInt::zero() {
            return self.inv().unwrap().pow(&(-exp));
        }
        let mut result = Self::one(&self.modulus);
        let mut base = self.clone();
        let mut e = exp.clone();
        while e > BigInt::zero() {
            if &e % BigInt::from(2) == BigInt::one() {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            e /= 2;
        }
        result
    }

    /// Convert to u64 if possible.
    pub fn to_u64(&self) -> Option<u64> {
        self.value.to_u64()
    }

    /// Get a reference to the inner BigInt value.
    pub fn as_bigint(&self) -> &BigInt {
        &self.value
    }
}

impl fmt::Display for FiniteField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_basic() {
        let a = FiniteField::from_u64(5);
        let b = FiniteField::from_u64(7);
        let c = a.add(&b);
        assert_eq!(c.to_u64(), Some(12));
    }

    #[test]
    fn sub_basic() {
        let a = FiniteField::from_u64(10);
        let b = FiniteField::from_u64(3);
        let c = a.sub(&b);
        assert_eq!(c.to_u64(), Some(7));
    }

    #[test]
    fn mul_basic() {
        let a = FiniteField::from_u64(6);
        let b = FiniteField::from_u64(7);
        let c = a.mul(&b);
        assert_eq!(c.to_u64(), Some(42));
    }

    #[test]
    fn inv_and_div() {
        let a = FiniteField::from_u64(3);
        let inv = a.inv().unwrap();
        let product = a.mul(&inv);
        assert_eq!(product.to_u64(), Some(1));
    }

    #[test]
    fn inv_zero_is_none() {
        let zero = FiniteField::from_u64(0);
        assert!(zero.inv().is_none());
    }

    #[test]
    fn pow_basic() {
        let a = FiniteField::from_u64(2);
        let result = a.pow(&BigInt::from(10));
        assert_eq!(result.to_u64(), Some(1024));
    }

    #[test]
    fn modular_reduction() {
        let modulus = BigInt::from(7u64);
        let a = FiniteField::new(BigInt::from(15), modulus);
        assert_eq!(a.to_u64(), Some(1)); // 15 % 7 = 1
    }

    #[test]
    fn negative_values() {
        let modulus = BigInt::from(7u64);
        let a = FiniteField::new(BigInt::from(-3), modulus);
        assert_eq!(a.to_u64(), Some(4)); // -3 mod 7 = 4
    }
}
