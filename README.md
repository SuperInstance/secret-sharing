# si-secret-sharing

Shamir's secret sharing, verifiable secret sharing (Feldman VSS), reconstruction
from subsets, and threshold variants over finite fields in Rust.

## Features

- **Finite field arithmetic** — modular arithmetic over large primes
- **Shamir's secret sharing** — (t, n) threshold scheme
- **Feldman VSS** — verifiable secret sharing with public commitments
- **Lagrange reconstruction** — recover secrets from any qualifying subset
- **Threshold utilities** — majority, any-two, and custom configurations

## Usage

```rust
use si_secret_sharing::{ShamirScheme, Reconstructor, FiniteField};

let secret = FiniteField::from_u64(42);
let shares = ShamirScheme::split(&secret, 5, 3); // 5 shares, threshold 3

let subset: Vec<&Share> = shares.iter().take(3).collect();
let recovered = Reconstructor::lagrange_interpolate(&subset);
assert_eq!(recovered, secret);
```

## License

MIT OR Apache-2.0
