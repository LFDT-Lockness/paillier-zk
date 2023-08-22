//! ZK-proof of Paillier-Blum modulus. Called Пmod or Rmod in the CGGMP21 paper.
//!
//! ## Description
//! A party P has a modulus `N = pq`, with p and q being Blum primes, and
//! `gcd(N, phi(N)) = 1`. P wants to prove that those equalities about N hold,
//! without disclosing p and q.
//!
//! ## Example
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use rug::{Integer, Complete};
//! let mut rng = rand_core::OsRng;
//!
//! // 0. Prover P derives two Blum primes and makes a Paillier-Blum modulus
//! let p = let p = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
//! let q = let p = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
//! let n = (&p * &q).complete();
//!
//! // 1. P computes a non-interactive proof that `n` is a Paillier-Blum modulus:
//! use paillier_zk::paillier_blum_modulus as p;
//!
//! // Security parameter
//! const SECURITY: usize = 33;
//! // Verifier and prover share the same state
//! let prover_shared_state = sha2::Sha256::default();
//! let verifier_shared_state = sha2::Sha256::default();
//!
//! let data = p::Data { n };
//! let pdata = p::PrivateData { p, q };
//!
//! let (commitment, proof) =
//!     p::non_interactive::prove::<{SECURITY}, _, _>(
//!         prover_shared_state,
//!         &data,
//!         &pdata,
//!         &mut rng,
//!     )?;
//!
//! // 2. P sends `data, commitment, proof` to the verifier V
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Proof<{SECURITY}>) { }
//! send(&data, &commitment, &proof);
//!
//! // 3. V receives and verifies the proof:
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//!
//! p::non_interactive::verify::<{SECURITY}, _>(
//!     verifier_shared_state,
//!     &data,
//!     &commitment,
//!     &proof,
//! )?;
//! # Ok(()) }
//! ```
//! If the verification succeeded, V can continue communication with P

use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Public data that both parties know: the Paillier-Blum modulus
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Data {
    pub n: Integer,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    pub p: Integer,
    pub q: Integer,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment {
    pub w: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
///
/// Consists of `M` singular challenges
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Challenge<const M: usize> {
    pub ys: [Integer; M],
}

/// A part of proof. Having enough of those guarantees security
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProofPoint {
    pub x: Integer,
    pub a: bool,
    pub b: bool,
    pub z: Integer,
}

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]. Consists of M proofs for each challenge
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof<const M: usize> {
    #[cfg_attr(
        // A trick to serialize arbitrary size arrays
        feature = "serde",
        serde(with = "serde_with::As::<[serde_with::Same; M]>")
    )]
    pub points: [ProofPoint; M],
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::common::sqrt::{blum_sqrt, find_residue, sample_non_residue_in};
    use crate::{Error, ErrorReason, IntegerExt, InvalidProof, InvalidProofReason};

    use super::{Challenge, Commitment, Data, PrivateData, Proof, ProofPoint};

    /// Create random commitment
    pub fn commit<R: RngCore>(Data { ref n }: &Data, rng: &mut R) -> Commitment {
        Commitment {
            w: sample_non_residue_in(n, rng),
        }
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<const M: usize>(
        Data { ref n }: &Data,
        PrivateData { ref p, ref q }: &PrivateData,
        Commitment { ref w }: &Commitment,
        challenge: &Challenge<M>,
    ) -> Result<Proof<M>, Error> {
        let sqrt = |x| blum_sqrt(&x, p, q, n);
        let phi = (p - 1u8).complete() * (q - 1u8).complete();
        let n_inverse = n.invert_ref(&phi).ok_or(ErrorReason::Invert)?.into();

        let points = challenge.ys.clone().map(|y| {
            let z = y.pow_mod_ref(&n_inverse, n)?.into();
            let (a, b, y_) = find_residue(&y, w, p, q, n);
            let x = sqrt(sqrt(y_));
            Some(ProofPoint { x, a, b, z })
        });
        if points.iter().any(Option::is_none) {
            return Err(ErrorReason::ModPow.into());
        }
        // we checked that all points are `Some(_)`. Have to do this trick as array::try_map is
        // not yet in stable
        let points = points.map(Option::unwrap);
        Ok(Proof { points })
    }

    /// Verify the proof. If this succeeds, the relation Rmod holds with chance
    /// `1/2^M`
    pub fn verify<const M: usize>(
        data: &Data,
        commitment: &Commitment,
        challenge: &Challenge<M>,
        proof: &Proof<M>,
    ) -> Result<(), InvalidProof> {
        if data.n.is_probably_prime(25) != rug::integer::IsPrime::No {
            return Err(InvalidProofReason::ModulusIsPrime.into());
        }
        if data.n.is_even() {
            return Err(InvalidProofReason::ModulusIsEven.into());
        }
        for (point, y) in proof.points.iter().zip(challenge.ys.iter()) {
            if Integer::from(
                point
                    .z
                    .pow_mod_ref(&data.n, &data.n)
                    .ok_or(InvalidProofReason::ModPow)?,
            ) != *y
            {
                return Err(InvalidProofReason::IncorrectNthRoot.into());
            }
            let y = y.clone();
            let y = if point.a { &data.n - y } else { y };
            let y = if point.b {
                (y * &commitment.w).modulo(&data.n)
            } else {
                y
            };
            if Integer::from(
                point
                    .x
                    .pow_mod_ref(&4.into(), &data.n)
                    .ok_or(InvalidProofReason::ModPow)?,
            ) != y
            {
                return Err(InvalidProofReason::IncorrectFourthRoot.into());
            }
        }
        Ok(())
    }

    /// Generate random challenge
    ///
    /// `data` parameter is used to generate challenge in correct range
    pub fn challenge<const M: usize, R: RngCore>(
        Data { ref n }: &Data,
        rng: &mut R,
    ) -> Challenge<M> {
        let ys = [(); M].map(|()| {
            n.random_below_ref(&mut fast_paillier::utils::external_rand(rng))
                .into()
        });
        Challenge { ys }
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use crate::{Error, InvalidProof};

    use super::{Challenge, Commitment, Data, PrivateData, Proof};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<const M: usize, R: RngCore, D>(
        shared_state: D,
        data: &Data,
        pdata: &PrivateData,
        rng: &mut R,
    ) -> Result<(Commitment, Proof<M>), Error>
    where
        D: Digest<OutputSize = U32> + Clone,
    {
        let commitment = super::interactive::commit(data, rng);
        let challenge = challenge(shared_state, data, &commitment);
        let proof = super::interactive::prove(data, pdata, &commitment, &challenge)?;
        Ok((commitment, proof))
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<const M: usize, D>(
        shared_state: D,
        data: &Data,
        commitment: &Commitment,
        proof: &Proof<M>,
    ) -> Result<(), InvalidProof>
    where
        D: Digest<OutputSize = U32> + Clone,
    {
        let challenge = challenge(shared_state, data, commitment);
        super::interactive::verify(data, commitment, &challenge, proof)
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<const M: usize, D>(
        shared_state: D,
        Data { ref n }: &Data,
        commitment: &Commitment,
    ) -> Challenge<M>
    where
        D: Digest,
    {
        let shared_state = shared_state.finalize();
        let hash = |d: D| {
            let order = rug::integer::Order::Msf;
            d.chain_update(&shared_state)
                .chain_update(n.to_digits::<u8>(order))
                .chain_update(commitment.w.to_digits::<u8>(order))
                .finalize()
        };
        let mut rng = crate::common::rng::HashRng::new(hash);
        // since we can't use Default and Integer isn't copy, we initialize
        // like this
        let ys = [(); M].map(|()| {
            n.random_below_ref(&mut fast_paillier::utils::external_rand(&mut rng))
                .into()
        });
        Challenge { ys }
    }
}

#[cfg(test)]
mod test {
    use rug::{Complete, Integer};

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let p = blum_prime(256, &mut rng);
        let q = blum_prime(256, &mut rng);
        let n = (&p * &q).complete();
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let shared_state = sha2::Sha256::default();
        let (commitment, proof) = super::non_interactive::prove::<65, _, _>(
            shared_state.clone(),
            &data,
            &pdata,
            &mut rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &data, &commitment, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{e:?}"),
        }
    }

    #[test]
    fn failing() {
        let mut rng = rand_dev::DevRng::new();
        let p = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
        let q = loop {
            // non blum prime
            let q = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
            if q.mod_u(4) == 1 {
                break q;
            }
        };
        let n = (&p * &q).complete();
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let shared_state = sha2::Sha256::default();
        let (commitment, proof) = super::non_interactive::prove::<65, _, _>(
            shared_state.clone(),
            &data,
            &pdata,
            &mut rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &data, &commitment, &proof);
        if r.is_ok() {
            panic!("proof should not pass");
        }
    }

    fn blum_prime<R: rand_core::RngCore>(s: u32, rng: &mut R) -> Integer {
        loop {
            let p = fast_paillier::utils::generate_safe_prime(rng, s);
            if p.mod_u(4) == 3 {
                break p;
            }
        }
    }
}
