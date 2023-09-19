//! ZK-proof for factoring of a RSA modulus. Called Пfac or Rfac in the CGGMP21
//! paper.
//!
//! ## Description
//!
//! A party P has a modulus `N = pq`. P wants to prove to a verifier V that p
//! and q are sufficiently large without disclosing p or q, with p and q each no
//! larger than `sqrt(N) * 2^l`, or equivalently no smaller than `sqrt(N) /
//! 2^l`
//!
//! ## Example
//!
//! ```rust
//! use rug::{Integer, Complete};
//! use paillier_zk::no_small_factor::non_interactive as p;
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         verifier_aux: p::Aux,
//! #     );
//! # }
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let shared_state_prover = sha2::Sha256::default();
//! let shared_state_verifier = sha2::Sha256::default();
//! let mut rng = rand_core::OsRng;
//! # let mut rng = rand_dev::DevRng::new();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters, and
//! // agree on the level of security
//!
//! let aux: p::Aux = pregenerated::verifier_aux();
//! let security = p::SecurityParams {
//!     l: 4,
//!     epsilon: 128,
//!     q: (Integer::ONE << 128_u32).complete(),
//! };
//!
//! // 1. Prover prepares the data to obtain proof about
//!
//! let p = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
//! let q = fast_paillier::utils::generate_safe_prime(&mut rng, 256);
//! let n = (&p * &q).complete();
//! let n_root = n.sqrt_ref().complete();
//! let data = p::Data {
//!     n: &n,
//!     n_root: &n_root,
//! };
//!
//! // 2. Prover computes a non-interactive proof that both factors are large enough
//!
//! let proof = p::prove(
//!     shared_state_prover,
//!     &aux,
//!     data,
//!     p::PrivateData { p: &p, q: &q },
//!     &security,
//!     rng,
//! )?;
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &Integer, _: &p::Proof) {  }
//! send(data.n, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data.n, proof);
//! let (n, proof) = recv();
//! let n_root = n.sqrt_ref().complete();;
//! let data = p::Data {
//!     n: &n,
//!     n_root: &n_root,
//! };
//! p::verify(shared_state_verifier, &aux, data, &security, &proof)?;
//! # Ok(()) }
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use rug::Integer;

pub use crate::common::{Aux, InvalidProof};

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, security parameter for bit size of plaintext: it needs to
    /// differ from sqrt(n) not more than by 2^l
    pub l: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy)]
pub struct Data<'a> {
    /// N0 - rsa modulus
    pub n: &'a Integer,
    /// A number close to square root of n
    pub n_root: &'a Integer,
}

/// Private data of prover
#[derive(Debug, Clone, Copy)]
pub struct PrivateData<'a> {
    pub p: &'a Integer,
    pub q: &'a Integer,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub beta: Integer,
    pub mu: Integer,
    pub nu: Integer,
    pub r: Integer,
    pub x: Integer,
    pub y: Integer,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment {
    pub p: Integer,
    pub q: Integer,
    pub a: Integer,
    pub b: Integer,
    pub t: Integer,
    pub sigma: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

/// The ZK proof, computed by [`interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: Integer,
    pub z2: Integer,
    pub w1: Integer,
    pub w2: Integer,
    pub v: Integer,
}

/// Interactive version of the proof
pub mod interactive {
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::{
        common::{fail_if, fail_if_ne, IntegerExt, InvalidProofReason},
        Error,
    };

    use super::{
        Aux, Challenge, Commitment, Data, InvalidProof, PrivateCommitment, PrivateData, Proof,
        SecurityParams,
    };

    /// Create random commitment
    pub fn commit<R: RngCore>(
        aux: &Aux,
        data: Data,
        pdata: PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment, PrivateCommitment), Error> {
        let two_to_l = (Integer::ONE << security.l).complete();
        let two_to_l_plus_e = (Integer::ONE << (security.l + security.epsilon)).complete();
        let n_root_modulo = (&two_to_l_plus_e * data.n_root).complete();
        let l_n_circ_modulo = (&two_to_l * &aux.rsa_modulo).complete();
        let l_e_n_circ_modulo = (&two_to_l_plus_e * &aux.rsa_modulo).complete();
        let n_n_circ = (&aux.rsa_modulo * data.n).complete();

        let alpha = Integer::from_rng_pm(&n_root_modulo, &mut rng);
        let beta = Integer::from_rng_pm(&n_root_modulo, &mut rng);
        let mu = Integer::from_rng_pm(&l_n_circ_modulo, &mut rng);
        let nu = Integer::from_rng_pm(&l_n_circ_modulo, &mut rng);
        let sigma = Integer::from_rng_pm(&(&two_to_l * &n_n_circ).complete(), &mut rng);
        let r = Integer::from_rng_pm(&(&two_to_l_plus_e * &n_n_circ).complete(), &mut rng);
        let x = Integer::from_rng_pm(&l_e_n_circ_modulo, &mut rng);
        let y = Integer::from_rng_pm(&l_e_n_circ_modulo, &mut rng);

        let p = aux.combine(pdata.p, &mu)?;
        let q = aux.combine(pdata.q, &nu)?;
        let a = aux.combine(&alpha, &x)?;
        let b = aux.combine(&beta, &y)?;
        let t = aux.rsa_modulo.combine(&q, &alpha, &aux.t, &r)?;

        let commitment = Commitment {
            p,
            q,
            a,
            b,
            t,
            sigma,
        };
        let private_commitment = PrivateCommitment {
            alpha,
            beta,
            mu,
            nu,
            r,
            x,
            y,
        };
        Ok((commitment, private_commitment))
    }

    /// Generate random challenge
    ///
    /// `security` parameter is used to generate challenge in correct range
    pub fn challenge<R: RngCore>(security: &SecurityParams, rng: &mut R) -> Challenge {
        Integer::from_rng_pm(&security.q, rng)
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        pdata: PrivateData,
        comm: &Commitment,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        let sigma_circ = (&comm.sigma - &pcomm.nu * pdata.p).complete();

        Ok(Proof {
            z1: (&pcomm.alpha + challenge * pdata.p).complete(),
            z2: (&pcomm.beta + challenge * pdata.q).complete(),
            w1: (&pcomm.x + challenge * &pcomm.mu).complete(),
            w2: (&pcomm.y + challenge * &pcomm.nu).complete(),
            v: &pcomm.r + challenge * sigma_circ,
        })
    }

    /// Verify the proof
    pub fn verify(
        aux: &Aux,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        // check 1
        {
            let lhs = aux.combine(&proof.z1, &proof.w1)?;
            let rhs =
                aux.rsa_modulo
                    .combine(&commitment.a, Integer::ONE, &commitment.p, challenge)?;
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        // check 2
        {
            let lhs = aux.combine(&proof.z2, &proof.w2)?;
            let rhs =
                aux.rsa_modulo
                    .combine(&commitment.b, Integer::ONE, &commitment.q, challenge)?;
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        // check 3
        {
            let r = aux.combine(data.n, &commitment.sigma)?;
            let lhs = aux
                .rsa_modulo
                .combine(&commitment.q, &proof.z1, &aux.t, &proof.v)?;
            let rhs = aux
                .rsa_modulo
                .combine(&commitment.t, Integer::ONE, &r, challenge)?;
            fail_if_ne(InvalidProofReason::EqualityCheck(3), lhs, rhs)?;
        }
        let range = (Integer::from(1) << (security.l + security.epsilon)) * data.n_root;
        // range check for z1
        fail_if(InvalidProofReason::RangeCheck(1), proof.z1.is_in_pm(&range))?;
        // range check for z2
        fail_if(InvalidProofReason::RangeCheck(2), proof.z2.is_in_pm(&range))?;

        Ok(())
    }
}

/// Non-interactive version of the proof
pub mod non_interactive {
    use digest::{typenum::U32, Digest};
    use rand_core::RngCore;

    pub use crate::{Error, InvalidProof};

    pub use super::{Aux, Challenge, Data, PrivateData, SecurityParams};

    /// The ZK proof, computed by [`prove`]
    #[derive(Debug, Clone)]
    #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
    pub struct Proof {
        commitment: super::Commitment,
        proof: super::Proof,
    }

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<R: RngCore, D>(
        shared_state: D,
        aux: &Aux,
        data: Data,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: R,
    ) -> Result<Proof, Error>
    where
        D: Digest<OutputSize = U32>,
    {
        let (commitment, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge(shared_state, aux, data, &commitment, security);
        let proof = super::interactive::prove(pdata, &commitment, &pcomm, &challenge)?;
        Ok(Proof { commitment, proof })
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<D>(
        shared_state: D,
        aux: &Aux,
        data: Data,
        commitment: &super::Commitment,
        security: &SecurityParams,
    ) -> Challenge
    where
        D: Digest,
    {
        let shared_state = shared_state.finalize();
        let hash = |d: D| {
            let order = rug::integer::Order::Msf;
            d.chain_update(&shared_state)
                .chain_update(aux.s.to_digits::<u8>(order))
                .chain_update(aux.t.to_digits::<u8>(order))
                .chain_update(aux.rsa_modulo.to_digits::<u8>(order))
                .chain_update(data.n.to_digits::<u8>(order))
                .chain_update(data.n_root.to_digits::<u8>(order))
                .chain_update(commitment.p.to_digits::<u8>(order))
                .chain_update(commitment.q.to_digits::<u8>(order))
                .chain_update(commitment.a.to_digits::<u8>(order))
                .chain_update(commitment.b.to_digits::<u8>(order))
                .chain_update(commitment.t.to_digits::<u8>(order))
                .chain_update(commitment.sigma.to_digits::<u8>(order))
                .finalize()
        };
        let mut rng = crate::common::rng::HashRng::new(hash);
        super::interactive::challenge(security, &mut rng)
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<D>(
        shared_state: D,
        aux: &Aux,
        data: Data,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof>
    where
        D: Digest<OutputSize = U32>,
    {
        let challenge = challenge(shared_state, aux, data, &proof.commitment, security);
        super::interactive::verify(
            aux,
            data,
            &proof.commitment,
            security,
            &challenge,
            &proof.proof,
        )
    }
}

#[cfg(test)]
mod test {
    use rug::{Complete, Integer};

    use crate::common::test::generate_blum_prime;
    use crate::common::InvalidProofReason;

    // If q > 2^epsilon, the proof will never pass. We can make l however small
    // we wish though, provided the statement we want to prove holds

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let p = generate_blum_prime(&mut rng, 256);
        let q = generate_blum_prime(&mut rng, 256);
        let n = (&p * &q).complete();
        let n_root = n.sqrt_ref().complete();
        let data = super::Data {
            n: &n,
            n_root: &n_root,
        };
        let security = super::SecurityParams {
            l: 64,
            epsilon: 128,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let aux = crate::common::test::aux(&mut rng);
        let shared_state = sha2::Sha256::default();
        let proof = super::non_interactive::prove(
            shared_state.clone(),
            &aux,
            data,
            super::PrivateData { p: &p, q: &q },
            &security,
            rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &aux, data, &security, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("Proof should not fail with {e:?}"),
        }
    }

    #[test]
    fn failing() {
        let mut rng = rand_dev::DevRng::new();
        let p = generate_blum_prime(&mut rng, 128);
        let q = generate_blum_prime(&mut rng, 384);
        let n = (&p * &q).complete();
        let n_root = n.sqrt_ref().complete();
        let data = super::Data {
            n: &n,
            n_root: &n_root,
        };
        let security = super::SecurityParams {
            l: 4,
            epsilon: 128,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let aux = crate::common::test::aux(&mut rng);
        let shared_state = sha2::Sha256::default();
        let proof = super::non_interactive::prove(
            shared_state.clone(),
            &aux,
            data,
            super::PrivateData { p: &p, q: &q },
            &security,
            rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &aux, data, &security, &proof)
            .expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::RangeCheck(2) => (),
            e => panic!("Proof should not fail with {e:?}"),
        }
    }

    #[test]
    fn test_sqrt() {
        assert_eq!(Integer::from(1).sqrt(), Integer::from(1));
        assert_eq!(Integer::from(2).sqrt(), Integer::from(1));
        assert_eq!(Integer::from(3).sqrt(), Integer::from(1));
        assert_eq!(Integer::from(4).sqrt(), Integer::from(2));
        assert_eq!(Integer::from(5).sqrt(), Integer::from(2));
        assert_eq!(Integer::from(6).sqrt(), Integer::from(2));
        assert_eq!(Integer::from(7).sqrt(), Integer::from(2));
        assert_eq!(Integer::from(8).sqrt(), Integer::from(2));
        assert_eq!(Integer::from(9).sqrt(), Integer::from(3));
    }
}
