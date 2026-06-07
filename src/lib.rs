//! # si-secret-sharing
//!
//! Shamir's secret sharing over finite fields, verifiable secret sharing
//! (Feldman VSS), reconstruction from subsets, and threshold variants.
//!
//! ## Modules
//!
//! - [`field`] — Finite field arithmetic (modular arithmetic over primes)
//! - [`shamir`] — Shamir's secret sharing scheme
//! - [`feldman`] — Feldman verifiable secret sharing
//! - [`reconstruct`] — Lagrange interpolation for secret reconstruction
//! - [`threshold`] — Threshold secret sharing utilities

pub mod field;
pub mod shamir;
pub mod feldman;
pub mod reconstruct;
pub mod threshold;

pub use field::FiniteField;
pub use shamir::ShamirScheme;
pub use feldman::FeldmanVSS;
pub use reconstruct::Reconstructor;
pub use threshold::ThresholdScheme;
