//! ZK-proof, called Пlog* or Rlog* in the CGGMP21 paper.
//!
//! ## Description
//!
//! A party P has a number `X = g ^ x`, with g being a generator of
//! multiplicative group G. P has encrypted x as C. P shares X and C with V and
//! wants to prove that the logarithm of X is the plaintext of C, and that the
//! plaintext (i.e. x) is at most L+1 bits.
//!
//! Given:
//! - `key0`, `pkey0` - pair of public and private keys in paillier cryptosystem
//! - `G` - a group of order `q` with generator `g`
//! - `X = g ^ x` and `C = key0.encrypt(x)` - data to obtain proof about
//!
//! Prove:
//! - `decrypt(C) = log X`
//! - `bitsize(x) <= L`
//!
//! Disclosing only: `key0`, `C`, `X`
//!
//! ## Example
//!
//! ```no_run
//! # use paillier_zk::unknown_order::BigNumber;
//! use paillier_zk::group_element_vs_paillier_encryption_in_range as p;
//! use generic_ec::hash_to_curve::Tag;
//!
//! // Prover and verifier have a shared protocol state
//! let shared_state_prover = sha2::Sha256::default();
//! let shared_state_verifier = sha2::Sha256::default();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
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
//! // 1. Setup: prover prepares the paillier keys
//!
//! let private_key = libpaillier::DecryptionKey::random().unwrap();
//! let key0 = libpaillier::EncryptionKey::from(&private_key);
//!
//! // 2. Setup: prover has some plaintext, encrypts it and computes X
//!
//! type C = generic_ec_curves::Secp256r1;
//! let g = generic_ec::Point::<C>::generator().into();
//!
//! let plaintext: BigNumber = 228.into();
//! let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
//! let power = g * paillier_zk::convert_scalar(&plaintext);
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let mut rng = rand_core::OsRng::default();
//!
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 128,
//!     q: BigNumber::prime_from_rng(128, &mut rng),
//! };
//!
//! let data = p::Data { key0, c: ciphertext, x: power, g };
//! let pdata = p::PrivateData { x: plaintext, nonce };
//! let (commitment, proof) =
//!     p::non_interactive::prove(shared_state_prover, &aux, &data, &pdata, &security, rng).expect("proof failed");
//!
//! // 4. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<C: Curve>(_: &p::Data<C>, _: &p::Commitment<C>, _: &p::Proof) { todo!() }
//! # fn recv<C: Curve>() -> (p::Data<C>, p::Commitment<C>, p::Proof) { todo!() }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! let (data, commitment, proof) = recv::<C>();
//! p::non_interactive::verify(shared_state_verifier, &aux, &data, &commitment, &security, &proof);
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use generic_ec::{Curve, Point};
use libpaillier::{Ciphertext, EncryptionKey, Nonce};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::unknown_order::BigNumber;

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, bit size of +-plaintext
    pub l: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: BigNumber,
}

/// Public data that both parties know
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Data<C: Curve> {
    /// N0 in paper, public key that C was encrypted on
    pub key0: EncryptionKey,
    /// C in paper, logarithm of X encrypted on N0
    pub c: Ciphertext,
    /// X in paper, exponent of plaintext of C
    pub x: Point<C>,
    /// A generator in group
    pub g: Point<C>,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    /// x in paper, logarithm of X and plaintext of C
    pub x: BigNumber,
    /// rho in paper, nonce in encryption x -> C
    pub nonce: Nonce,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Commitment<C: Curve> {
    pub s: BigNumber,
    pub a: Ciphertext,
    pub y: Point<C>,
    pub d: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub mu: BigNumber,
    pub r: Nonce,
    pub gamma: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = BigNumber;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub z3: BigNumber,
}

pub use crate::common::Aux;

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use generic_ec::Curve;
    use libpaillier::unknown_order::BigNumber;
    use rand_core::RngCore;

    use crate::common::{combine, convert_scalar, gen_inversible, InvalidProof, ProtocolError};

    use super::{
        Aux, Challenge, Commitment, Data, PrivateCommitment, PrivateData, Proof, SecurityParams,
    };

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        aux: &Aux,
        data: &Data<C>,
        pdata: &PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment<C>, PrivateCommitment), ProtocolError> {
        let two_to_l = BigNumber::one() << (security.l + 1);
        let two_to_l_e = BigNumber::one() << (security.l + security.epsilon + 1);
        let modulo_l = two_to_l * &aux.rsa_modulo;
        let modulo_l_e = &two_to_l_e * &aux.rsa_modulo;

        let alpha = BigNumber::from_rng(&two_to_l_e, &mut rng);
        let mu = BigNumber::from_rng(&modulo_l, &mut rng);
        let r = gen_inversible(data.key0.n(), &mut rng);
        let gamma = BigNumber::from_rng(&modulo_l_e, &mut rng);

        let (a, _) = data
            .key0
            .encrypt(alpha.to_bytes(), Some(r.clone()))
            .ok_or(ProtocolError::EncryptionFailed)?;

        let commitment = Commitment {
            s: combine(&aux.s, &pdata.x, &aux.t, &mu, &aux.rsa_modulo),
            a,
            y: data.g * convert_scalar(&alpha),
            d: combine(&aux.s, &alpha, &aux.t, &gamma, &aux.rsa_modulo),
        };
        let private_commitment = PrivateCommitment {
            alpha,
            mu,
            r,
            gamma,
        };
        Ok((commitment, private_commitment))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<C: Curve>(
        data: &Data<C>,
        pdata: &PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Proof {
        Proof {
            z1: &pcomm.alpha + challenge * &pdata.x,
            z2: combine(
                &pcomm.r,
                &BigNumber::one(),
                &pdata.nonce,
                challenge,
                data.key0.n(),
            ),
            z3: &pcomm.gamma + challenge * &pcomm.mu,
        }
    }

    /// Verify the proof
    pub fn verify<C: Curve>(
        aux: &Aux,
        data: &Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        let one = BigNumber::one();
        fn fail_if(b: bool, msg: InvalidProof) -> Result<(), InvalidProof> {
            if b {
                Ok(())
            } else {
                Err(msg)
            }
        }
        // Three equality checks and one range check
        {
            let (lhs, _) = data
                .key0
                .encrypt(proof.z1.to_bytes(), Some(proof.z2.clone()))
                .ok_or(InvalidProof::EncryptionFailed)?;
            let rhs = combine(&commitment.a, &one, &data.c, challenge, data.key0.nn());
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(1))?;
        }
        {
            let lhs = data.g * convert_scalar(&proof.z1);
            let rhs = commitment.y + data.x * convert_scalar(challenge);
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(2))?;
        }
        {
            let lhs = combine(&aux.s, &proof.z1, &aux.t, &proof.z3, &aux.rsa_modulo);
            let rhs = combine(
                &commitment.d,
                &one,
                &commitment.s,
                challenge,
                &aux.rsa_modulo,
            );
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(3))?;
        }
        fail_if(
            proof.z1 <= one << (security.l + security.epsilon + 1),
            InvalidProof::RangeCheckFailed(4),
        )?;

        Ok(())
    }

    /// Generate random challenge
    ///
    /// `data` parameter is used to generate challenge in correct range
    pub fn challenge<R>(rng: &mut R, security: &SecurityParams) -> BigNumber
    where
        R: RngCore,
    {
        // double the range to account for +-
        let m = BigNumber::from(2) * &security.q;
        BigNumber::from_rng(&m, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use generic_ec::{hash_to_curve::FromHash, Curve, Scalar};
    use libpaillier::unknown_order::BigNumber;
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use crate::common::{InvalidProof, ProtocolError};

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<C: Curve, R: RngCore, D: Digest>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        pdata: &PrivateData,
        security: &SecurityParams,
        rng: R,
    ) -> Result<(Commitment<C>, Proof), ProtocolError>
    where
        Scalar<C>: generic_ec::hash_to_curve::FromHash,
        D: Digest<OutputSize = U32>,
    {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge);
        Ok((comm, proof))
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<C: Curve, D: Digest>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof>
    where
        Scalar<C>: generic_ec::hash_to_curve::FromHash,
        D: Digest<OutputSize = U32>,
    {
        let challenge = challenge(shared_state, aux, data, commitment, security);
        super::interactive::verify(aux, data, commitment, security, &challenge, proof)
    }

    /// Internal function for deriving challenge from protocol values
    /// deterministically
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
    ) -> Challenge
    where
        Scalar<C>: FromHash,
        D: Digest<OutputSize = U32>,
    {
        use rand_core::SeedableRng;
        let seed = shared_state
            .chain_update(aux.s.to_bytes())
            .chain_update(aux.t.to_bytes())
            .chain_update(aux.rsa_modulo.to_bytes())
            .chain_update(data.key0.to_bytes())
            .chain_update(data.c.to_bytes())
            .chain_update(data.x.to_bytes(true))
            .chain_update(commitment.s.to_bytes())
            .chain_update(commitment.a.to_bytes())
            .chain_update(commitment.y.to_bytes(true))
            .chain_update(commitment.d.to_bytes())
            .chain_update((security.l as u64).to_le_bytes())
            .chain_update((security.epsilon as u64).to_le_bytes())
            .finalize();

        let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed.into());
        let m = BigNumber::from(2) * &security.q;
        BigNumber::from_rng(&m, &mut rng)
    }
}

#[cfg(test)]
mod test {
    use generic_ec::{hash_to_curve::FromHash, Curve, Scalar};
    use libpaillier::unknown_order::BigNumber;

    use crate::common::test::{nonce, random_key};
    use crate::common::{convert_scalar, InvalidProof};

    fn run<R: rand_core::RngCore, C: Curve>(
        mut rng: R,
        security: super::SecurityParams,
        plaintext: BigNumber,
    ) -> Result<(), crate::common::InvalidProof>
    where
        Scalar<C>: FromHash,
    {
        let private_key0 = random_key(&mut rng).unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let nonce = nonce(&mut rng, key0.n());
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), nonce).unwrap();
        let g = generic_ec::Point::<C>::generator() * generic_ec::Scalar::<C>::from(1337);
        let x = g * convert_scalar(&plaintext);

        let data = super::Data {
            key0,
            c: ciphertext,
            x,
            g,
        };
        let pdata = super::PrivateData {
            x: plaintext,
            nonce,
        };

        let aux = crate::common::test::aux(&mut rng);

        let shared_state = sha2::Sha256::default();

        let (commitment, proof) = super::non_interactive::prove(
            shared_state.clone(),
            &aux,
            &data,
            &pdata,
            &security,
            rng,
        )
        .unwrap();

        super::non_interactive::verify(shared_state, &aux, &data, &commitment, &security, &proof)
    }

    fn passing_test<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_core::OsRng::default();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 256,
            q: BigNumber::prime_from_rng(128, &mut rng),
        };
        let plaintext = BigNumber::from(228);
        let r = run(rng, security, plaintext);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{e:?}"),
        }
    }

    fn failing_test<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_core::OsRng::default();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
            q: BigNumber::prime_from_rng(128, &mut rng),
        };
        let plaintext = BigNumber::from(1) << (security.l + security.epsilon + 1);
        let r = run(rng, security, plaintext);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(InvalidProof::RangeCheckFailed(_)) => (),
            Err(e) => panic!("proof should not fail with: {e:?}"),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec_curves::rust_crypto::Secp256r1>()
    }
    #[test]
    fn failing_p256() {
        failing_test::<generic_ec_curves::rust_crypto::Secp256r1>()
    }

    #[test]
    fn passing_million() {
        passing_test::<crate::curve::C>()
    }
    #[test]
    fn failing_million() {
        failing_test::<crate::curve::C>()
    }

    // see notes in
    // [crate::paillier_encryption_in_range::test::rejected_with_probability_1_over_2]
    // for motivation and design of the following test.
    // Altough no security estimate was given in the paper, my own calculations
    // show that the parameters here achieve the probability about as good as in
    // other tests

    #[test]
    fn mul_rejected_with_probability_1_over_2() {
        use rand_core::SeedableRng;
        fn maybe_rejected(mut rng: rand_chacha::ChaCha20Rng) -> bool {
            let security = super::SecurityParams {
                l: 1024,
                epsilon: 130,
                q: BigNumber::prime_from_rng(128, &mut rng),
            };
            let plaintext = (BigNumber::from(1) << (security.l + 1)) - 1;
            let r = run::<_, generic_ec_curves::rust_crypto::Secp256r1>(rng, security, plaintext);
            match r {
                Ok(()) => true,
                Err(crate::common::InvalidProof::RangeCheckFailed(4)) => false,
                Err(e) => panic!("proof should not fail with: {e:?}"),
            }
        }

        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(2);
        assert!(!maybe_rejected(rng), "should fail");
        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(3);
        assert!(maybe_rejected(rng), "should pass");
    }
}
