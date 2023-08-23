#![deny(clippy::disallowed_methods)]

use thiserror::Error;

mod common;
pub mod group_element_vs_paillier_encryption_in_range;
pub mod no_small_factor;
pub mod paillier_affine_operation_in_range;
pub mod paillier_blum_modulus;
pub mod paillier_encryption_in_range;

#[cfg(test)]
mod curve;

use common::InvalidProofReason;
pub use common::{BadExponent, IntegerExt, InvalidProof, PaillierError};

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
