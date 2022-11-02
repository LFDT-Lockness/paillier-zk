pub mod sqrt;

use crate::unknown_order::BigNumber;

/// Generate element in Zm*. Does so by trial.
pub fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, mut rng: R) -> BigNumber {
    loop {
        let r = BigNumber::from_rng(modulo, &mut rng);
        if r.gcd(modulo) == 1.into() {
            break r;
        }
    }
}

/// Compute l^le * r^re modulo m
pub fn combine(
    l: &BigNumber,
    le: &BigNumber,
    r: &BigNumber,
    re: &BigNumber,
    m: &BigNumber,
) -> BigNumber {
    l.modpow(le, m).modmul(&r.modpow(re, m), m)
}

/// Reason for failure. Mainly interesting for debugging purposes
#[derive(Debug, PartialEq, Eq)]
pub enum InvalidProof {
    /// One equality doesn't hold. Parameterized by equality index
    EqualityCheckFailed(usize),
    /// One range check doesn't hold. Parameterized by check index
    RangeCheckFailed(usize),
    /// Encryption of supplied data failed when attempting to verify
    EncryptionFailed,
}

#[derive(Debug, PartialEq, Eq)]
pub enum ProtocolError {
    /// Encryption of supplied data failed when computing proof
    EncryptionFailed,
    /// Hashing of supplied data failed when computing proof or challenge
    HashFailed,
}
