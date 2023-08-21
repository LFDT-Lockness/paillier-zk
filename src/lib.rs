#![deny(clippy::disallowed_methods)]

// Don't forget to add new backends as they appear in unknown_order
#[cfg(all(not(feature = "gmp"), not(feature = "rust")))]
compile_error!(
    "Attempting to build paillier-zk without a bignumber backend specified: \
    please enable `rust` or `gmp` feature of `paillier-zk` (`openssl` backend is \
    disabled as it doesn't support deterministic RNG)"
);

use thiserror::Error;

mod common;
pub mod group_element_vs_paillier_encryption_in_range;
pub mod no_small_factor;
pub mod paillier_affine_operation_in_range;
pub mod paillier_blum_modulus;
pub mod paillier_encryption_in_range;
// pub mod paillier_decryption_modulo_q;

#[cfg(test)]
mod curve;

/// Underlying paillier library for which the proofs are made. Use this to get
/// the correct version of the library
pub use libpaillier;
/// Underlying big number implementation. Use this to get
/// the correct version of the library
pub use libpaillier::unknown_order;

use common::InvalidProofReason;
pub use common::{
    BadExponent, BigNumberExt, IntegerExt, InvalidProof, PaillierError, SafePaillierDecryptionExt,
    SafePaillierEncryptionExt,
};

/// Library general error type
#[derive(Debug, Error)]
#[error(transparent)]
pub struct Error(#[from] ErrorReason);

#[derive(Debug, Error)]
enum ErrorReason {
    #[error("couldn't evaluate modpow")]
    ModPow,
    #[error("couldn't encrypt a message")]
    Encryption,
    #[error("can't find multiplicative inverse")]
    Invert,
    #[error("paillier error")]
    Paillier(#[source] fast_paillier::Error),
}

impl From<BadExponent> for Error {
    fn from(_err: BadExponent) -> Self {
        Error(ErrorReason::ModPow)
    }
}

impl From<PaillierError> for Error {
    fn from(_err: PaillierError) -> Self {
        Error(ErrorReason::Encryption)
    }
}

impl From<fast_paillier::Error> for Error {
    fn from(err: fast_paillier::Error) -> Self {
        Self(ErrorReason::Paillier(err))
    }
}
