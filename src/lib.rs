mod common;
pub mod group_element_vs_paillier_encryption_in_range;
pub mod paillier_affine_operation_in_range;
pub mod paillier_blum_modulus;
pub mod paillier_encryption_in_range;
pub mod paillier_decryption_modulo_q;

pub(crate) mod curve;

/// Underlying paillier library for which the proofs are made. Use this to get
/// the correct version of the library
pub use libpaillier;
/// Underlying big number implementation. Use this to get
/// the correct version of the library
pub use libpaillier::unknown_order;

/// Bit size in Пenc and Пaff-g
/// TODO: choose appropriate value
pub const L: usize = 228;
/// Bit size overshoot in Пenc and Пaff-g
/// TODO: choose appropriate value
pub const EPSILON: usize = 322;
/// Challenges amount in Пmod
/// TODO: choose appropriate value
pub const M: usize = 13;
