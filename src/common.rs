pub mod sqrt;

use crate::unknown_order::BigNumber;

/// Auxiliary data known to both prover and verifier
pub struct Aux {
    /// ring-pedersen parameter
    pub s: BigNumber,
    /// ring-pedersen parameter
    pub t: BigNumber,
    /// N^ in paper
    pub rsa_modulo: BigNumber,
}

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
pub(crate) fn combine(
    l: &BigNumber,
    le: &BigNumber,
    r: &BigNumber,
    re: &BigNumber,
    m: &BigNumber,
) -> BigNumber {
    l.modpow(le, m).modmul(&r.modpow(re, m), m)
}

/// Embed BigInt into chosen scalar type
pub fn convert_scalar<C: generic_ec::Curve>(x: &BigNumber) -> generic_ec::Scalar<C> {
    generic_ec::Scalar::<C>::from_be_bytes_mod_order(x.to_bytes())
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

#[cfg(test)]
mod test {
    use libpaillier::unknown_order::BigNumber;

    #[test]
    fn conversion() {
        // checks that bignumbers use BE encoding, same as the method we use in
        // conversion
        type Scalar = generic_ec::Scalar<generic_ec_curves::Secp256r1>;
        let number: u64 = 0x11_22_33_44_55_66_77_88;
        let bignumber = BigNumber::from(number);
        let scalar1 = Scalar::from(number);
        let scalar2 = super::convert_scalar(&bignumber);
        assert_eq!(scalar1, scalar2);
    }
}
