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

/// Reason for failure. If the proof failes, you should only be interested in a
/// reason for debugging purposes
#[derive(Debug, PartialEq, Eq)]
pub enum InvalidProof {
    /// One equality doesn't hold. Parameterized by equality index
    EqualityCheckFailed(usize),
    /// One range check doesn't hold. Parameterized by check index
    RangeCheckFailed(usize),
    /// Encryption of supplied data failed when attempting to verify
    EncryptionFailed,
}

/// Unexpeted error that can happen in a protocol. You should probably panic if
/// you see this.
#[derive(Debug, PartialEq, Eq)]
pub enum ProtocolError {
    /// Encryption of supplied data failed when computing proof
    EncryptionFailed,
    /// Hashing of supplied data failed when computing proof or challenge
    HashFailed,
}

/// Regular paillier encryption methods are easy to misuse and generate
/// an undeterministic nonce, we replace them with those functions
pub trait SafePaillierExt {
    fn encrypt_with<M: AsRef<[u8]>>(
        &self,
        x: M,
        nonce: libpaillier::Nonce,
    ) -> Option<libpaillier::Ciphertext>;
    fn encrypt_with_random<M: AsRef<[u8]>, R: rand_core::RngCore>(
        &self,
        x: M,
        rng: &mut R,
    ) -> Option<(libpaillier::Ciphertext, libpaillier::Nonce)>;
}

impl SafePaillierExt for libpaillier::EncryptionKey {
    fn encrypt_with<M: AsRef<[u8]>>(
        &self,
        x: M,
        nonce: libpaillier::Nonce,
    ) -> Option<libpaillier::Ciphertext> {
        #[allow(clippy::disallowed_methods)]
        self.encrypt(x, Some(nonce)).map(|(e, _)| e)
    }

    fn encrypt_with_random<M: AsRef<[u8]>, R: rand_core::RngCore>(
        &self,
        x: M,
        rng: &mut R,
    ) -> Option<(libpaillier::Ciphertext, libpaillier::Nonce)> {
        let nonce = libpaillier::Nonce::from_rng(self.n(), rng);
        #[allow(clippy::disallowed_methods)]
        self.encrypt(x, Some(nonce))
    }
}

pub trait BigNumberExt {
    /// Generate element in Zm*. Does so by trial.
    fn gen_inversible<R: rand_core::RngCore>(rng: &mut R, modulo: &BigNumber) -> Self;

    /// Compute l^le * r^re modulo m
    fn combine(l: &Self, le: &Self, r: &Self, re: &Self, m: &Self) -> Self;

    /// Embed BigInt into chosen scalar type
    fn to_scalar<C: generic_ec::Curve>(&self) -> generic_ec::Scalar<C>;

    /// Generates a random integer in interval `[-range; range]`
    fn from_rng_pm<R: rand_core::RngCore>(rng: &mut R, range: &Self) -> Self;
}

impl BigNumberExt for BigNumber {
    fn gen_inversible<R: rand_core::RngCore>(rng: &mut R, modulo: &BigNumber) -> Self {
        loop {
            let r = BigNumber::from_rng(modulo, rng);
            if r.gcd(modulo) == 1.into() {
                break r;
            }
        }
    }

    fn combine(l: &Self, le: &Self, r: &Self, re: &Self, m: &Self) -> Self {
        l.modpow(le, m).modmul(&r.modpow(re, m), m)
    }

    fn to_scalar<C: generic_ec::Curve>(&self) -> generic_ec::Scalar<C> {
        generic_ec::Scalar::<C>::from_be_bytes_mod_order(self.to_bytes())
    }

    fn from_rng_pm<R: rand_core::RngCore>(rng: &mut R, range: &Self) -> Self {
        let n = BigNumber::from_rng(&(range * 2), rng);
        n - range
    }
}

#[cfg(test)]
pub mod test {
    use libpaillier::unknown_order::BigNumber;

    use super::BigNumberExt;

    #[test]
    fn conversion() {
        // checks that bignumbers use BE encoding, same as the method we use in
        // conversion
        type Scalar = generic_ec::Scalar<generic_ec_curves::Secp256r1>;
        let number: u64 = 0x11_22_33_44_55_66_77_88;
        let bignumber = BigNumber::from(number);
        let scalar1 = Scalar::from(number);
        let scalar2 = bignumber.to_scalar();
        assert_eq!(scalar1, scalar2);
    }

    pub fn random_key<R: rand_core::RngCore>(rng: &mut R) -> Option<libpaillier::DecryptionKey> {
        let p = BigNumber::prime_from_rng(1024, rng);
        let q = BigNumber::prime_from_rng(1024, rng);
        libpaillier::DecryptionKey::with_primes_unchecked(&p, &q)
    }

    pub fn aux<R: rand_core::RngCore>(rng: &mut R) -> super::Aux {
        let p = BigNumber::prime_from_rng(1024, rng);
        let q = BigNumber::prime_from_rng(1024, rng);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        super::Aux { s, t, rsa_modulo }
    }
}
