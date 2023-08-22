pub mod rng;
pub mod sqrt;

use generic_ec::Scalar;
use rug::{Complete, Integer};

use crate::unknown_order::BigNumber;

pub mod for_rug {
    use rug::Integer;

    /// Auxiliary data known to both prover and verifier
    pub struct Aux {
        /// ring-pedersen parameter
        pub s: Integer,
        /// ring-pedersen parameter
        pub t: Integer,
        /// N^ in paper
        pub rsa_modulo: Integer,
    }
}

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
    #[error("paillier encryption failed")]
    PaillierEnc,
    #[error("paillier homomorphic op failed")]
    PaillierOp,
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

impl From<PaillierError> for InvalidProof {
    fn from(_err: PaillierError) -> Self {
        InvalidProof(InvalidProofReason::Encryption)
    }
}

/// Regular paillier encryption methods are easy to misuse and generate
/// an undeterministic nonce, we replace them with those functions
pub trait SafePaillierEncryptionExt {
    /// Encrypts an integer `x` in `[-N/2; N/2)` with `nonce` in `Z*_N`
    ///
    /// Returns error if inputs are not in specified range
    fn encrypt_with(
        &self,
        x: &BigNumber,
        nonce: &libpaillier::Nonce,
    ) -> Result<libpaillier::Ciphertext, PaillierError>;

    /// Encrypts an integer `x` in `[-N/2; N/2)`, generates a random nonce
    /// from provided PRNG
    ///
    /// Returns error if `x` is not in the specified range
    fn encrypt_with_random<R: rand_core::RngCore>(
        &self,
        x: &BigNumber,
        rng: &mut R,
    ) -> Result<(libpaillier::Ciphertext, libpaillier::Nonce), PaillierError>;

    /// Homomorphic addition of two ciphertexts
    fn oadd(
        &self,
        c1: &libpaillier::Ciphertext,
        c2: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError>;

    /// Homomorphic multiplication of scalar (should be unsigned integer) at ciphertext
    fn omul(
        &self,
        scalar: &BigNumber,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError>;

    /// Homomorphic negation of a ciphertext
    fn oneg(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError>;
}

pub trait SafePaillierDecryptionExt {
    /// Decrypts ciphertext to plaintext in `[-N/2; N/2)`
    ///
    /// Returns error if ciphertext is malformed
    fn decrypt_to_bigint(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<BigNumber, PaillierError>;
}

impl SafePaillierEncryptionExt for libpaillier::EncryptionKey {
    #[allow(clippy::disallowed_methods)]
    fn encrypt_with(
        &self,
        x: &BigNumber,
        nonce: &libpaillier::Nonce,
    ) -> Result<libpaillier::Ciphertext, PaillierError> {
        let x_twice = x << 1;
        if !(-self.n() <= x_twice && &x_twice < self.n()) {
            return Err(PaillierError);
        }

        let x = x.nmod(self.n());
        if nonce.gcd(self.n()) != 1.into() {
            return Err(PaillierError);
        }

        // a = (1 + N)^x mod N^2 = (1 + xN) mod N^2
        let a = BigNumber::one().modadd(&(&x * self.n()), self.nn());
        // b = nonce^N mod N^2
        let b = nonce.modpow(self.n(), self.nn());

        let c = a.modmul(&b, self.nn());
        Ok(c)
    }

    fn encrypt_with_random<R: rand_core::RngCore>(
        &self,
        x: &BigNumber,
        rng: &mut R,
    ) -> Result<(libpaillier::Ciphertext, libpaillier::Nonce), PaillierError> {
        let nonce = libpaillier::Nonce::from_rng(self.n(), rng);
        self.encrypt_with(x, &nonce).map(|ciph| (ciph, nonce))
    }

    fn oadd(
        &self,
        c1: &libpaillier::Ciphertext,
        c2: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError> {
        #[allow(clippy::disallowed_methods)]
        self.add(c1, c2).ok_or(PaillierError)
    }

    fn omul(
        &self,
        scalar: &BigNumber,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError> {
        #[allow(clippy::disallowed_methods)]
        if scalar >= &BigNumber::zero() {
            self.mul(ciphertext, scalar).ok_or(PaillierError)
        } else {
            self.mul(&self.oneg(ciphertext)?, &-scalar)
                .ok_or(PaillierError)
        }
    }

    fn oneg(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<libpaillier::Ciphertext, PaillierError> {
        ciphertext.invert(self.nn()).ok_or(PaillierError)
    }
}

impl SafePaillierDecryptionExt for libpaillier::DecryptionKey {
    fn decrypt_to_bigint(
        &self,
        ciphertext: &libpaillier::Ciphertext,
    ) -> Result<BigNumber, PaillierError> {
        #![allow(clippy::disallowed_methods)]
        let plaintext = self.decrypt(ciphertext).ok_or(PaillierError)?;
        let plaintext = BigNumber::from_slice(plaintext);
        if &(&plaintext << 1) >= self.n() {
            Ok(plaintext - self.n())
        } else {
            Ok(plaintext)
        }
    }
}

/// Error indicating that encryption failed
#[derive(Clone, Copy, Debug, thiserror::Error)]
#[error("paillier encryption failed")]
pub struct PaillierError;

pub trait BigNumberExt: Sized {
    /// Generate element in Zm*. Does so by trial.
    fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, rng: &mut R) -> Self;

    /// Compute l^le * r^re modulo self
    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent>;

    /// Embed BigInt into chosen scalar type
    fn to_scalar<C: generic_ec::Curve>(&self) -> Scalar<C>;

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

    fn to_scalar<C: generic_ec::Curve>(&self) -> Scalar<C> {
        if self >= &BigNumber::zero() {
            Scalar::<C>::from_be_bytes_mod_order(self.to_bytes())
        } else {
            -Scalar::<C>::from_be_bytes_mod_order((-self).to_bytes())
        }
    }

    fn curve_order<C: generic_ec::Curve>() -> BigNumber {
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

pub trait IntegerExt: Sized {
    /// Generate element in Zm*. Does so by trial.
    fn gen_inversible<R: rand_core::RngCore>(modulo: &Self, rng: &mut R) -> Self;

    /// Compute l^le * r^re modulo self
    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent>;

    /// Embed BigInt into chosen scalar type
    fn to_scalar<C: generic_ec::Curve>(&self) -> Scalar<C>;

    /// Returns prime order of curve C
    fn curve_order<C: generic_ec::Curve>() -> Self;

    /// Generates a random integer in interval `[-range; range]`
    fn from_rng_pm<R: rand_core::RngCore>(range: &Self, rng: &mut R) -> Self;

    /// Checks whether `self` is in interval `[-range; range]`
    fn is_in_pm(&self, range: &Self) -> bool;

    /// Returns `self mod n`
    fn modulo(&self, n: &Self) -> Self;

    /// Returns `self smod n` in `{-n/2, .., n/2}`
    fn signed_modulo(&self, n: &Self) -> Self;
}

impl IntegerExt for Integer {
    fn gen_inversible<R: rand_core::RngCore>(modulo: &Integer, rng: &mut R) -> Self {
        fast_paillier::utils::sample_in_mult_group(rng, modulo)
    }

    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent> {
        let l_to_le: Integer = l.pow_mod_ref(le, self).ok_or(BadExponent)?.into();
        let r_to_re: Integer = r.pow_mod_ref(re, self).ok_or(BadExponent)?.into();
        Ok((l_to_le * r_to_re).modulo(self))
    }

    fn to_scalar<C: generic_ec::Curve>(&self) -> Scalar<C> {
        let bytes_be = self.to_digits::<u8>(rug::integer::Order::Msf);
        let s = Scalar::<C>::from_be_bytes_mod_order(bytes_be);
        if self.cmp0().is_ge() {
            s
        } else {
            -s
        }
    }

    fn curve_order<C: generic_ec::Curve>() -> Self {
        let order_minus_one = -Scalar::<C>::one();
        let i = Integer::from_digits(&order_minus_one.to_be_bytes(), rug::integer::Order::Msf);
        i + 1
    }

    fn from_rng_pm<R: rand_core::RngCore>(range: &Self, rng: &mut R) -> Self {
        let mut rng = fast_paillier::utils::external_rand(rng);
        let range_twice = range.clone() << 1u32;
        range_twice.random_below(&mut rng) - range
    }

    fn is_in_pm(&self, range: &Self) -> bool {
        let minus_range = -range.clone();
        minus_range <= *self && self <= range
    }

    fn modulo(&self, n: &Self) -> Self {
        fast_paillier::utils::IntegerExt::modulo(self, n)
    }

    fn signed_modulo(&self, n: &Self) -> Self {
        let self_mod_n = self.modulo(n);
        let half_n = (n >> 1_u32).complete();
        if half_n.is_odd() && self_mod_n <= half_n || self_mod_n < half_n {
            self_mod_n
        } else {
            self_mod_n - n
        }
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

/// A common logic shared across tests and doctests
#[cfg(any(test, feature = "__internal_doctest"))]
pub mod test {
    use libpaillier::unknown_order::BigNumber;
    use rug::{Complete, Integer};

    use super::IntegerExt;

    pub fn random_key<R: rand_core::RngCore>(rng: &mut R) -> Option<libpaillier::DecryptionKey> {
        let p = BigNumber::prime_from_rng(1024, rng);
        let q = BigNumber::prime_from_rng(1024, rng);
        libpaillier::DecryptionKey::with_primes_unchecked(&p, &q)
    }

    pub fn random_key_rug<R: rand_core::RngCore>(
        rng: &mut R,
    ) -> Option<fast_paillier::DecryptionKey> {
        let p = generate_blum_prime(rng, 1024);
        let q = generate_blum_prime(rng, 1024);
        fast_paillier::DecryptionKey::from_primes(p, q).ok()
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

    pub fn aux_rug<R: rand_core::RngCore>(rng: &mut R) -> super::for_rug::Aux {
        let p = generate_blum_prime(rng, 1024);
        let q = generate_blum_prime(rng, 1024);
        let n = (&p * &q).complete();

        let (s, t) = {
            let phi_n = (p.clone() - 1u8) * (q.clone() - 1u8);
            let r = Integer::gen_inversible(&n, rng);
            let lambda = phi_n.random_below(&mut fast_paillier::utils::external_rand(rng));

            let t = r.square().modulo(&n);
            let s = t.pow_mod_ref(&lambda, &n).unwrap().into();

            (s, t)
        };

        super::for_rug::Aux {
            s,
            t,
            rsa_modulo: n,
        }
    }

    pub fn generate_blum_prime(rng: &mut impl rand_core::RngCore, bits_size: u32) -> Integer {
        loop {
            let n = generate_prime(rng, bits_size);
            if n.mod_u(4) == 3 {
                break n;
            }
        }
    }

    pub fn generate_prime(rng: &mut impl rand_core::RngCore, bits_size: u32) -> Integer {
        let mut n: Integer =
            Integer::random_bits(bits_size, &mut fast_paillier::utils::external_rand(rng)).into();
        n.set_bit(bits_size - 1, true);
        n.next_prime_mut();
        n
    }
}

#[cfg(test)]
mod _test {
    use std::iter;

    use libpaillier::unknown_order::BigNumber;
    use rand_dev::DevRng;
    use rug::{Complete, Integer};

    use crate::SafePaillierEncryptionExt;

    use super::{test::random_key, BigNumberExt, IntegerExt, SafePaillierDecryptionExt};

    #[test]
    fn pailler_encryption_decryption() {
        let mut rng = DevRng::new();
        let dk = random_key(&mut rng).unwrap();
        let ek = libpaillier::EncryptionKey::from(&dk);

        let corner_cases = (-10..=10).map(|i| (ek.n() / 2) + i);
        let random_xs = iter::repeat_with(|| BigNumber::from_rng(ek.n(), &mut rng)).take(20);
        let xs = corner_cases.chain(random_xs).collect::<Vec<_>>();

        for x in xs {
            let x = if &x * 2 >= *ek.n() { x - ek.n() } else { x };

            let (ciphertext, _) = ek.encrypt_with_random(&x, &mut rng).unwrap();
            let decrypted = dk.decrypt_to_bigint(&ciphertext).unwrap();

            assert_eq!(x, decrypted);
        }
    }

    #[test]
    fn conversion() {
        // checks that bignumbers use BE encoding, same as the method we use in
        // conversion
        type Scalar = generic_ec::Scalar<generic_ec::curves::Secp256r1>;
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

    #[test]
    fn to_scalar_encoding() {
        type E = generic_ec::curves::Secp256k1;

        let bytes = [123u8, 231u8];
        let int = u16::from_be_bytes(bytes);
        let bn = rug::Integer::from(int);
        let scalar = bn.to_scalar();
        assert_eq!(scalar, generic_ec::Scalar::<E>::from(int));

        assert_eq!(bn.to_digits::<u8>(rug::integer::Order::Msf), &bytes);

        let curve_order = Integer::curve_order::<E>();
        assert_eq!(curve_order.to_scalar(), generic_ec::Scalar::<E>::zero());
        assert_eq!(
            (curve_order - 1u8).to_scalar(),
            -generic_ec::Scalar::<E>::one()
        );
    }

    #[test]
    fn signed_modulo() {
        let n = Integer::from(7);

        assert_eq!(Integer::from(0).signed_modulo(&n), 0);
        assert_eq!(Integer::from(1).signed_modulo(&n), 1);
        assert_eq!(Integer::from(2).signed_modulo(&n), 2);
        assert_eq!(Integer::from(3).signed_modulo(&n), 3);
        assert_eq!(Integer::from(4).signed_modulo(&n), -3);
        assert_eq!(Integer::from(5).signed_modulo(&n), -2);
        assert_eq!(Integer::from(6).signed_modulo(&n), -1);
        assert_eq!(Integer::from(7).signed_modulo(&n), 0);
        assert_eq!(Integer::from(8).signed_modulo(&n), 1);

        let n = Integer::from(4);
        assert_eq!(Integer::from(0).signed_modulo(&n), 0);
        assert_eq!(Integer::from(1).signed_modulo(&n), 1);
        assert_eq!(Integer::from(2).signed_modulo(&n), -2);
        assert_eq!(Integer::from(3).signed_modulo(&n), -1);
    }
}
