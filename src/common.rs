pub mod rng;
pub mod sqrt;

use std::sync::Arc;

use generic_ec::Scalar;
use rug::{Complete, Integer};

/// Auxiliary data known to both prover and verifier
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug)]
pub struct Aux {
    /// ring-pedersen parameter
    pub s: Integer,
    /// ring-pedersen parameter
    pub t: Integer,
    /// N^ in paper
    pub rsa_modulo: Integer,
    /// Precomuted table for computing `s^x t^y mod rsa_modulo` faster
    ///
    /// If absent, optimization is disabled.
    ///
    /// We wrap the table into [`Arc`] to make cloning cheaper. Note that during serialization/deserialization
    /// process, serde may create some extra copies of the table.
    pub multiexp: Option<Arc<crate::multiexp::MultiexpTable>>,
}

impl Aux {
    pub fn combine(&self, x: &Integer, y: &Integer) -> Result<Integer, BadExponent> {
        if let Some(table) = &self.multiexp {
            match table.prod_exp(x, y) {
                Some(res) => return Ok(res),
                None if cfg!(debug_assertions) => {
                    return Err(BadExponentReason::PrecompTable.into())
                }
                None => {
                    // When debug assertions are disabled, we fallback to naive exponentiation
                }
            }
        }

        self.combine_naive(x, y)
    }

    fn combine_naive(&self, x: &Integer, y: &Integer) -> Result<Integer, BadExponent> {
        let s_to_x: Integer = self
            .s
            .pow_mod_ref(x, &self.rsa_modulo)
            .ok_or(BadExponentReason::Undefined)?
            .into();
        let t_to_y: Integer = self
            .t
            .pow_mod_ref(y, &self.rsa_modulo)
            .ok_or(BadExponentReason::Undefined)?
            .into();
        Ok((s_to_x * t_to_y).modulo(&self.rsa_modulo))
    }
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

/// Error indicating that encryption failed
#[derive(Clone, Copy, Debug, thiserror::Error)]
#[error("paillier encryption failed")]
pub struct PaillierError;

pub trait IntegerExt: Sized {
    /// Generate element in Zm*. Does so by trial.
    fn gen_invertible<R: rand_core::RngCore>(modulo: &Self, rng: &mut R) -> Self;

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

    /// Returns `self smod n`
    ///
    /// For odd `n`, result is in `{-n/2, .., n/2}`. For even `n`, result is in
    /// `{-n/2, .., n/2 - 1}`
    fn signed_modulo(&self, n: &Self) -> Self;
}

impl IntegerExt for Integer {
    fn gen_invertible<R: rand_core::RngCore>(modulo: &Integer, rng: &mut R) -> Self {
        fast_paillier::utils::sample_in_mult_group(rng, modulo)
    }

    fn combine(&self, l: &Self, le: &Self, r: &Self, re: &Self) -> Result<Self, BadExponent> {
        let l_to_le: Integer = l
            .pow_mod_ref(le, self)
            .ok_or(BadExponentReason::Undefined)?
            .into();
        let r_to_re: Integer = r
            .pow_mod_ref(re, self)
            .ok_or(BadExponentReason::Undefined)?
            .into();
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

    fn signed_modulo(&self, n: &Self) -> Self {
        let self_mod_n = self.modulo_ref(n).complete();
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
#[error(transparent)]
pub struct BadExponent(#[from] BadExponentReason);

impl BadExponent {
    /// Constructs an error that exponent is undefined
    pub fn undefined() -> Self {
        Self(BadExponentReason::Undefined)
    }
}

#[derive(Clone, Copy, Debug, thiserror::Error)]
enum BadExponentReason {
    #[error("exponent is undefined")]
    Undefined,
    #[error("multiexp error: precomputation table could not perform multiexponentiation")]
    PrecompTable,
}

/// Returns `Err(err)` if `assertion` is false
pub fn fail_if<E>(err: E, assertion: bool) -> Result<(), E> {
    if assertion {
        Ok(())
    } else {
        Err(err)
    }
}

/// Returns `Err(err)` if `lhs != rhs`
pub fn fail_if_ne<T: PartialEq, E>(err: E, lhs: T, rhs: T) -> Result<(), E> {
    if lhs == rhs {
        Ok(())
    } else {
        Err(err)
    }
}

/// A common logic shared across tests and doctests
#[cfg(test)]
pub mod test {
    use rug::{Complete, Integer};

    use super::IntegerExt;

    pub fn random_key<R: rand_core::RngCore>(rng: &mut R) -> Option<fast_paillier::DecryptionKey> {
        let p = generate_blum_prime(rng, 1024);
        let q = generate_blum_prime(rng, 1024);
        fast_paillier::DecryptionKey::from_primes(p, q).ok()
    }

    pub fn aux<R: rand_core::RngCore>(rng: &mut R) -> super::Aux {
        let p = generate_blum_prime(rng, 1024);
        let q = generate_blum_prime(rng, 1024);
        let n = (&p * &q).complete();

        let (s, t) = {
            let phi_n = (p.clone() - 1u8) * (q.clone() - 1u8);
            let r = Integer::gen_invertible(&n, rng);
            let lambda = phi_n.random_below(&mut fast_paillier::utils::external_rand(rng));

            let t = r.square().modulo(&n);
            let s = t.pow_mod_ref(&lambda, &n).unwrap().into();

            (s, t)
        };

        super::Aux {
            s,
            t,
            rsa_modulo: n,
            multiexp: None,
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
    use rug::Integer;

    use super::IntegerExt;

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
