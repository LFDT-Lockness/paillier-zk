//! ZK-proof for factoring of a RSA modulus. Called Пfac or Rfac in the CGGMP21
//! paper.
//!
//! ## Description
//!
//! A party P has a modulus `N = pq`. P wants to prove to a verifier V that p
//! and q are sufficiently large without disclosing p or q, with p and q each no
//! larger than `sqrt(N) * 2^(l+1)`, or equivalently no smaller than `sqrt(N) /
//! 2^(l+1)`
//!
//! ## Example
//!
//! ```no_run
//! # use paillier_zk::unknown_order::BigNumber;
//! # fn sqrt(x: &BigNumber) -> BigNumber { todo!() }
//! use paillier_zk::no_small_factor::non_interactive as p;
//! let shared_state_prover = sha2::Sha256::default();
//! let shared_state_verifier = sha2::Sha256::default();
//! let mut rng = rand_core::OsRng::default();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters, and
//! // agree on the level of security
//!
//! let p = BigNumber::prime(1024);
//! let q = BigNumber::prime(1024);
//! let rsa_modulo = p * q;
//! let s: BigNumber = 123.into();
//! let t: BigNumber = 321.into();
//! assert_eq!(s.gcd(&rsa_modulo), 1.into());
//! assert_eq!(t.gcd(&rsa_modulo), 1.into());
//!
//! let aux = p::Aux { s, t, rsa_modulo };
//!
//! let security = p::SecurityParams {
//!     l: 4,
//!     epsilon: 128,
//!     q: BigNumber::prime_from_rng(128, &mut rng),
//! };
//!
//! // 1. Prover prepares the data to obtain proof about
//!
//! let p = BigNumber::prime_from_rng(256, &mut rng);
//! let q = BigNumber::prime_from_rng(256, &mut rng);
//! let n = &p * &q;
//! let n_root = sqrt(&n);
//! let data = p::Data {
//!     n: &n,
//!     n_root: &n_root,
//! };
//!
//! // 2. Prover computes a non-interactive proof that both factors are small
//!
//! let proof = p::prove(
//!     shared_state_prover,
//!     &aux,
//!     data,
//!     p::PrivateData { p: &p, q: &q },
//!     &security,
//!     rng,
//! );
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &BigNumber, _: &p::Proof) { todo!() }
//! # fn recv() -> (BigNumber, p::Proof) { todo!() }
//! send(data.n, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! let (n, proof) = recv();
//! let n_root = sqrt(&n);
//! let data = p::Data {
//!     n: &n,
//!     n_root: &n_root,
//! };
//! p::verify(shared_state_verifier, &aux, data, &security, &proof);
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::unknown_order::BigNumber;

pub use crate::common::InvalidProof;

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, security parameter for bit size of plaintext: it needs to
    /// differ from sqrt(n) not more than by 2^(l+1)
    pub l: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: BigNumber,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy)]
pub struct Data<'a> {
    /// N0 - rsa modulus
    pub n: &'a BigNumber,
    /// A number close to square root of n
    pub n_root: &'a BigNumber,
}

/// Private data of prover
#[derive(Debug, Clone, Copy)]
pub struct PrivateData<'a> {
    pub p: &'a BigNumber,
    pub q: &'a BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub beta: BigNumber,
    pub mu: BigNumber,
    pub nu: BigNumber,
    pub r: BigNumber,
    pub x: BigNumber,
    pub y: BigNumber,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment {
    pub p: BigNumber,
    pub q: BigNumber,
    pub a: BigNumber,
    pub b: BigNumber,
    pub t: BigNumber,
    pub sigma: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = BigNumber;

/// The ZK proof, computed by [`interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub w1: BigNumber,
    pub w2: BigNumber,
    pub v: BigNumber,
}

pub use crate::common::Aux;

pub mod interactive {
    use crate::{common::combine, unknown_order::BigNumber};
    use rand_core::RngCore;

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
    ) -> (Commitment, PrivateCommitment) {
        // add 1 to exponents to account for +-
        let two_to_l = BigNumber::from(1) << (security.l + 1);
        let two_to_l_plus_e = BigNumber::from(1) << (security.l + security.epsilon + 1);
        let n_root_modulo = &two_to_l_plus_e * data.n_root;
        let l_n_circ_modulo = &two_to_l * &aux.rsa_modulo;
        let l_e_n_circ_modulo = &two_to_l_plus_e * &aux.rsa_modulo;
        let n_n_circ = &aux.rsa_modulo * data.n;

        let alpha = BigNumber::from_rng(&n_root_modulo, &mut rng);
        let beta = BigNumber::from_rng(&n_root_modulo, &mut rng);
        let mu = BigNumber::from_rng(&l_n_circ_modulo, &mut rng);
        let nu = BigNumber::from_rng(&l_n_circ_modulo, &mut rng);
        let sigma = BigNumber::from_rng(&(&two_to_l * &n_n_circ), &mut rng);
        let r = BigNumber::from_rng(&(&two_to_l_plus_e * &n_n_circ), &mut rng);
        let x = BigNumber::from_rng(&l_e_n_circ_modulo, &mut rng);
        let y = BigNumber::from_rng(&l_e_n_circ_modulo, &mut rng);

        let p = combine(&aux.s, pdata.p, &aux.t, &mu, &aux.rsa_modulo);
        let q = combine(&aux.s, pdata.q, &aux.t, &nu, &aux.rsa_modulo);
        let a = combine(&aux.s, &alpha, &aux.t, &x, &aux.rsa_modulo);
        let b = combine(&aux.s, &beta, &aux.t, &y, &aux.rsa_modulo);
        let t = combine(&q, &alpha, &aux.t, &r, &aux.rsa_modulo);

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
        (commitment, private_commitment)
    }

    /// Generate random challenge
    ///
    /// `security` parameter is used to generate challenge in correct range
    pub fn challenge<R: RngCore>(security: &SecurityParams, rng: &mut R) -> Challenge {
        // double the range to account for +-
        let m = BigNumber::from(2) * &security.q;
        BigNumber::from_rng(&m, rng)
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        pdata: PrivateData,
        comm: &Commitment,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Proof {
        let sigma_circ = &comm.sigma - &pcomm.nu * pdata.p;

        Proof {
            z1: &pcomm.alpha + challenge * pdata.p,
            z2: &pcomm.beta + challenge * pdata.q,
            w1: &pcomm.x + challenge * &pcomm.mu,
            w2: &pcomm.y + challenge * &pcomm.nu,
            v: &pcomm.r + challenge * sigma_circ,
        }
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
        let one = BigNumber::one();
        // check 1
        {
            let lhs = combine(&aux.s, &proof.z1, &aux.t, &proof.w1, &aux.rsa_modulo);
            let rhs = combine(
                &commitment.a,
                &one,
                &commitment.p,
                challenge,
                &aux.rsa_modulo,
            );
            if lhs != rhs {
                return Err(InvalidProof::EqualityCheckFailed(1));
            }
        }
        // check 2
        {
            let lhs = combine(&aux.s, &proof.z2, &aux.t, &proof.w2, &aux.rsa_modulo);
            let rhs = combine(
                &commitment.b,
                &one,
                &commitment.q,
                challenge,
                &aux.rsa_modulo,
            );
            if lhs != rhs {
                return Err(InvalidProof::EqualityCheckFailed(2));
            }
        }
        // check 3
        {
            let r = combine(&aux.s, data.n, &aux.t, &commitment.sigma, &aux.rsa_modulo);
            let lhs = combine(&commitment.q, &proof.z1, &aux.t, &proof.v, &aux.rsa_modulo);
            let rhs = combine(&commitment.t, &one, &r, challenge, &aux.rsa_modulo);
            if lhs != rhs {
                return Err(InvalidProof::EqualityCheckFailed(3));
            }
        }
        let range = (BigNumber::from(1) << (security.l + security.epsilon + 1)) * data.n_root;
        // range check for z1
        if proof.z1 > range {
            return Err(InvalidProof::RangeCheckFailed(1));
        }
        // range check for z2
        if proof.z2 > range {
            return Err(InvalidProof::RangeCheckFailed(2));
        }

        Ok(())
    }
}

pub mod non_interactive {
    use crate::unknown_order::BigNumber;
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    pub use crate::common::{InvalidProof, ProtocolError};

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
    ) -> Proof
    where
        D: Digest<OutputSize = U32>,
    {
        let (commitment, pcomm) = super::interactive::commit(aux, data, pdata, security, rng);
        let challenge = challenge(shared_state, aux, data, &commitment, security);
        let proof = super::interactive::prove(pdata, &commitment, &pcomm, &challenge);
        Proof { commitment, proof }
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
        D: Digest<OutputSize = U32>,
    {
        use rand_core::SeedableRng;
        let seed = shared_state
            .chain_update(aux.s.to_bytes())
            .chain_update(aux.t.to_bytes())
            .chain_update(aux.rsa_modulo.to_bytes())
            .chain_update(data.n.to_bytes())
            .chain_update(data.n_root.to_bytes())
            .chain_update(commitment.p.to_bytes())
            .chain_update(commitment.q.to_bytes())
            .chain_update(commitment.a.to_bytes())
            .chain_update(commitment.b.to_bytes())
            .chain_update(commitment.t.to_bytes())
            .chain_update(commitment.sigma.to_bytes())
            .finalize();
        let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed.into());
        let m = BigNumber::from(2) * &security.q;
        BigNumber::from_rng(&m, &mut rng)
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
    use crate::common::InvalidProof;
    use crate::unknown_order::BigNumber;

    // If q > 2^epsilon, the proof will never pass. We can make l however small
    // we wish though, provided the statement we want to prove holds

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let p = BigNumber::prime_from_rng(256, &mut rng);
        let q = BigNumber::prime_from_rng(256, &mut rng);
        let n = &p * &q;
        let n_root = sqrt(&n);
        let data = super::Data {
            n: &n,
            n_root: &n_root,
        };
        let security = super::SecurityParams {
            l: 4,
            epsilon: 128,
            q: BigNumber::prime_from_rng(128, &mut rng),
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
        );
        let r = super::non_interactive::verify(shared_state, &aux, data, &security, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("Proof should not fail with {e:?}"),
        }
    }

    #[test]
    fn failing() {
        let mut rng = rand_dev::DevRng::new();
        let p = BigNumber::prime_from_rng(128, &mut rng);
        let q = BigNumber::prime_from_rng(384, &mut rng);
        let n = &p * &q;
        let n_root = sqrt(&n);
        let data = super::Data {
            n: &n,
            n_root: &n_root,
        };
        let security = super::SecurityParams {
            l: 4,
            epsilon: 128,
            q: BigNumber::prime_from_rng(128, &mut rng),
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
        );
        let r = super::non_interactive::verify(shared_state, &aux, data, &security, &proof);
        match r {
            Ok(()) => panic!("Proof should not pass"),
            Err(InvalidProof::RangeCheckFailed(2)) => (),
            Err(e) => panic!("Proof should not fail with {e:?}"),
        }
    }

    #[test]
    fn test_sqrt() {
        assert_eq!(sqrt(&BigNumber::from(1)), BigNumber::from(1));
        assert_eq!(sqrt(&BigNumber::from(2)), BigNumber::from(1));
        assert_eq!(sqrt(&BigNumber::from(3)), BigNumber::from(1));
        assert_eq!(sqrt(&BigNumber::from(4)), BigNumber::from(2));
        assert_eq!(sqrt(&BigNumber::from(5)), BigNumber::from(2));
        assert_eq!(sqrt(&BigNumber::from(6)), BigNumber::from(2));
        assert_eq!(sqrt(&BigNumber::from(7)), BigNumber::from(2));
        assert_eq!(sqrt(&BigNumber::from(8)), BigNumber::from(2));
        assert_eq!(sqrt(&BigNumber::from(9)), BigNumber::from(3));
    }

    /// Binary search for square root
    fn sqrt(x: &BigNumber) -> BigNumber {
        let mut low = BigNumber::one();
        let mut high = x.clone();
        let mut count = 1000;
        while low < &high - 1 {
            let mid = (&high + &low) / 2;
            let test: BigNumber = &mid * &mid;
            match test.cmp(x) {
                std::cmp::Ordering::Equal => return mid,
                std::cmp::Ordering::Less => {
                    low = mid;
                }
                std::cmp::Ordering::Greater => {
                    high = mid;
                }
            }
            count -= 1;
            assert_ne!(count, 0);
        }
        low
    }
}
