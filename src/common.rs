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

/// Error indicating that proof is invalid
#[derive(Debug, Clone, thiserror::Error)]
#[error("invalid proof")]
pub struct InvalidProof(
    #[source]
    #[from]
    InvalidProofReason,
);

/// Reason for failure. If the proof failes, you should only be interested in a
/// reason for debugging purposes
#[derive(Debug, PartialEq, Eq, Clone, Copy, thiserror::Error)]
pub enum InvalidProofReason {
    /// One equality doesn't hold. Parameterized by equality index
    #[error("equality check failed {0}")]
    EqualityCheck(usize),
    /// One range check doesn't hold. Parameterized by check index
    #[error("range check failed {0}")]
    RangeCheck(usize),
    /// Encryption of supplied data failed when attempting to verify
    #[error("encryption failed")]
    Encryption,
    /// Failed to evaluate powmod
    #[error("powmod failed")]
    ModPow,
    /// Paillier-Blum modulus is prime
    #[error("modulus is prime")]
    ModulusIsPrime,
    /// Paillier-Blum modulus is even
    #[error("modulus is even")]
    ModulusIsEven,
    /// Proof's z value in n-th power does not equal commitment value
    #[error("incorrect nth root")]
    IncorrectNthRoot,
    /// Proof's x value in 4-th power does not equal commitment value
    #[error("incorrect 4th root")]
    IncorrectFourthRoot,
}

impl InvalidProof {
    #[cfg(test)]
    pub(crate) fn reason(&self) -> InvalidProofReason {
        self.0
    }
}

impl From<BadExponent> for InvalidProof {
    fn from(_err: BadExponent) -> Self {
        InvalidProofReason::ModPow.into()
    }
}

impl From<EncryptionError> for InvalidProof {
    fn from(_err: EncryptionError) -> Self {
        InvalidProof(InvalidProofReason::Encryption)
    }
}

/// Regular paillier encryption methods are easy to misuse and generate
/// an undeterministic nonce, we replace them with those functions
pub trait SafePaillierExt {
    fn encrypt_with(
        &self,
        x: &BigNumber,
        nonce: &libpaillier::Nonce,
    ) -> Result<libpaillier::Ciphertext, EncryptionError>;
    fn encrypt_with_random<R: rand_core::RngCore>(
        &self,
        x: &BigNumber,
        rng: &mut R,
    ) -> Result<(libpaillier::Ciphertext, libpaillier::Nonce), EncryptionError>;

    /// Homomorphic addition of two ciphertexts
    fn oadd(
        &self,
        c1: &libpaillier::Ciphertext,
        c2: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError>;

    /// Homomorphic multiplication of scalar (should be unsigned integer) at ciphertext
    fn omul(
        &self,
        scalar: &BigNumber,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError>;

    /// Homomorphic negation of a ciphertext
    fn oneg(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError>;
}

impl SafePaillierExt for libpaillier::EncryptionKey {
    fn encrypt_with(
        &self,
        x: &BigNumber,
        nonce: &libpaillier::Nonce,
    ) -> Result<libpaillier::Ciphertext, EncryptionError> {
        if !(&BigNumber::zero() <= x && x < self.n()) {
            return Err(EncryptionError);
        }
        #[allow(clippy::disallowed_methods)]
        self.encrypt(x.to_bytes(), Some(nonce.clone()))
            .map(|(e, _)| e)
            .ok_or(EncryptionError)
    }

    fn encrypt_with_random<R: rand_core::RngCore>(
        &self,
        x: &BigNumber,
        rng: &mut R,
    ) -> Result<(libpaillier::Ciphertext, libpaillier::Nonce), EncryptionError> {
        if !(&BigNumber::zero() <= x && x < self.n()) {
            return Err(EncryptionError);
        }
        let nonce = libpaillier::Nonce::from_rng(self.n(), rng);
        #[allow(clippy::disallowed_methods)]
        self.encrypt(x.to_bytes(), Some(nonce))
            .ok_or(EncryptionError)
    }

    fn oadd(
        &self,
        c1: &libpaillier::Ciphertext,
        c2: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError> {
        #[allow(clippy::disallowed_methods)]
        self.add(c1, c2).ok_or(EncryptionError)
    }

    fn omul(
        &self,
        scalar: &BigNumber,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError> {
        #[allow(clippy::disallowed_methods)]
        if scalar >= &BigNumber::zero() {
            self.mul(ciphertext, scalar).ok_or(EncryptionError)
        } else {
            self.mul(&self.oneg(ciphertext)?, &-scalar)
                .ok_or(EncryptionError)
        }
    }

    fn oneg(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, EncryptionError> {
        ciphertext.invert(self.nn()).ok_or(EncryptionError)
    }
}

/// Error indicating that encryption failed
///
/// Returned by [SafePaillierExt] methods
#[derive(Clone, Copy, Debug, thiserror::Error)]
#[error("paillier encryption failed")]
pub struct EncryptionError;

pub trait BigNumberExt: Sized {
    /// Generate element in Zm*. Does so by trial.
    fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, rng: &mut R) -> Self;

    /// Compute l^le * r^re modulo self
    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent>;

    /// Embed BigInt into chosen scalar type
    fn to_scalar<C: generic_ec::Curve>(&self) -> generic_ec::Scalar<C>;

    /// Returns prime order of curve C
    fn curve_order<C: generic_ec::Curve>() -> BigNumber;

    /// Generates a random integer in interval `[-range; range]`
    fn from_rng_pm<R: rand_core::RngCore>(range: &Self, rng: &mut R) -> Self;

    /// Computes self^exponent mod modulo
    ///
    /// Unlike [`BigNumber::modpow`], this method correctly handles negative exponent. `Err(_)`
    /// is returned if modpow cannot be computed.
    fn powmod(&self, exponent: &Self, modulo: &Self) -> Result<Self, BadExponent>;

    /// Checks whether `self` is in interval `[-range; range]`
    fn is_in_pm(&self, range: &Self) -> bool;
}

impl BigNumberExt for BigNumber {
    fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, rng: &mut R) -> Self {
        loop {
            let r = BigNumber::from_rng(modulo, rng);
            if r.gcd(modulo) == 1.into() {
                break r;
            }
        }
    }

    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent> {
        Ok(l.powmod(le, self)?.modmul(&r.powmod(re, self)?, self))
    }

    fn to_scalar<C: generic_ec::Curve>(&self) -> generic_ec::Scalar<C> {
        if self >= &BigNumber::zero() {
            generic_ec::Scalar::<C>::from_be_bytes_mod_order(self.to_bytes())
        } else {
            -generic_ec::Scalar::<C>::from_be_bytes_mod_order((-self).to_bytes())
        }
    }

    fn curve_order<C: generic_ec::Curve>() -> BigNumber {
        use generic_ec::Scalar;
        let n_minus_one = Scalar::<C>::zero() - Scalar::one();
        let bn = BigNumber::from_slice(n_minus_one.to_be_bytes());
        bn + 1
    }

    fn from_rng_pm<R: rand_core::RngCore>(range: &Self, rng: &mut R) -> Self {
        let n = BigNumber::from_rng(&(range * 2), rng);
        n - range
    }

    fn powmod(&self, exponent: &Self, modulo: &Self) -> Result<Self, BadExponent> {
        if modulo <= &BigNumber::zero() {
            return Err(BadExponent);
        }

        #[allow(clippy::disallowed_methods)]
        if exponent < &BigNumber::zero() {
            Ok(BigNumber::modpow(
                &self.invert(modulo).ok_or(BadExponent)?,
                &(-exponent),
                modulo,
            ))
        } else {
            Ok(BigNumber::modpow(self, exponent, modulo))
        }
    }

    fn is_in_pm(&self, range: &Self) -> bool {
        let neg_range = -range;

        let low = &neg_range <= self;
        let high = self <= range;

        low & high
    }
}

/// Error indicating that computation cannot be evaluated because of bad exponent
///
/// Returned by [`BigNumberExt::powmod`] and other functions that do exponentiation internally
#[derive(Clone, Copy, Debug, thiserror::Error)]
#[error("exponent is undefined")]
pub struct BadExponent;

/// Returns `Err(err)` if `assertion` is false
pub(crate) fn fail_if<E>(err: E, assertion: bool) -> Result<(), E> {
    if assertion {
        Ok(())
    } else {
        Err(err)
    }
}

/// Returns `Err(err)` if `lhs != rhs`
pub(crate) fn fail_if_ne<T: PartialEq, E>(err: E, lhs: T, rhs: T) -> Result<(), E> {
    if lhs == rhs {
        Ok(())
    } else {
        Err(err)
    }
}

#[cfg(test)]
pub mod test {
    use libpaillier::unknown_order::BigNumber;
    use rand_dev::DevRng;

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

    #[test]
    fn generates_negative_numbers() {
        let mut rng = DevRng::new();

        let n = BigNumber::from(100);
        let mut generated_neg_numbers = 0;

        for _ in 0..32 {
            let i = BigNumber::from_rng_pm(&n, &mut rng);
            if i < BigNumber::zero() {
                generated_neg_numbers += 1
            }
        }

        assert!(generated_neg_numbers > 0);
        println!("generated {generated_neg_numbers} negative numbers");
    }

    #[test]
    fn arithmetic_with_negative_numbers_works_fine() {
        let two = BigNumber::from(2);
        let neg_2 = BigNumber::zero() - &two;

        assert_eq!(two + &neg_2, BigNumber::zero());
        assert_eq!(BigNumber::from(3) + &neg_2, BigNumber::one());
        assert_eq!(
            BigNumber::from(5) * &neg_2,
            BigNumber::zero() - BigNumber::from(10)
        );

        assert_eq!(
            BigNumber::from(5).modmul(&neg_2, &BigNumber::from(17)),
            BigNumber::from(7)
        );
        assert_eq!(
            BigNumber::one().modadd(&neg_2, &BigNumber::from(17)),
            BigNumber::from(16)
        );
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
