//! ZK-proof of paillier operation with group commitment in range. Called Пaff-g
//! or Raff-g in the CGGMP21 paper.
//!
//! ## Description
//!
//! A party P performs a paillier affine operation with C, Y, and X
//! obtaining `D = C*X + Y`. `X` and `Y` are encrypted values of `x` and `y`. P
//! then wants to prove that `y` and `x` are at most `L` and `L'` bits,
//! correspondingly, and P doesn't want to disclose none of the plaintexts
//!
//! Given:
//! - `key0`, `pkey0`, `key1`, `pkey1` - pairs of public and private keys in
//!   paillier cryptosystem
//! - `nonce_y`, `nonce` - nonces in paillier encryption
//! - `x`, `y` - some numbers
//! - `q`, `g` such that `<g> = Zq*` - prime order group
//! - `C` is some ciphertext encrypted by `key0`
//! - `Y = key1.encrypt(y, nonce_y)`
//! - `X = g * x`
//! - `D = oadd(enc(y, nonce), omul(x, C))` where `enc`, `oadd` and `omul` are
//!   paillier encryption, homomorphic addition and multiplication with `key0`
//!
//! Prove:
//! - `bitsize(abs(x)) <= l_x`
//! - `bitsize(abs(y)) <= l_y`
//!
//! Disclosing only: `key0`, `key1`, `C`, `D`, `Y`, `X`
//!
//! ## Example
//!
//! ```rust
//! use paillier_zk::{paillier_affine_operation_in_range as p, IntegerExt};
//! use rug::{Integer, Complete};
//! use generic_ec::{Point, curves::Secp256k1 as E};
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         verifier_aux: p::Aux,
//! #         someone_encryption_key0: fast_paillier::EncryptionKey,
//! #         someone_encryption_key1: fast_paillier::EncryptionKey,
//! #     );
//! # }
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Prover and verifier have a shared protocol state
//! let shared_state_prover = sha2::Sha256::default();
//! let shared_state_verifier = sha2::Sha256::default();
//!
//! let mut rng = rand_core::OsRng;
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
//!
//! let aux: p::Aux = pregenerated::verifier_aux();
//! let security = p::SecurityParams {
//!     l_x: 256,
//!     l_y: 848,
//!     epsilon: 230,
//!     q: (Integer::ONE << 128_u32).complete(),
//! };
//!
//! // 1. Setup: prover prepares the paillier keys
//!
//! // C and D are encrypted by this key
//! let key0: fast_paillier::EncryptionKey = pregenerated::someone_encryption_key0();
//! // Y is encrypted using this key
//! let key1: fast_paillier::EncryptionKey = pregenerated::someone_encryption_key1();
//!
//! // C is some number encrypted using key0. Neither of parties
//! // need to know the plaintext
//! let ciphertext_c = Integer::gen_inversible(&key0.nn(), &mut rng);
//!
//! // 2. Setup: prover prepares all plaintexts
//!
//! // x in paper
//! let plaintext_x = Integer::from_rng_pm(
//!     &(Integer::ONE << security.l_x).complete(),
//!     &mut rng,
//! );
//! // y in paper
//! let plaintext_y = Integer::from_rng_pm(
//!     &(Integer::ONE << security.l_y).complete(),
//!     &mut rng,
//! );
//!
//! // 3. Setup: prover encrypts everything on correct keys and remembers some nonces
//!
//! // X in paper
//! let ciphertext_x = Point::<E>::generator() * plaintext_x.to_scalar();
//! // Y and ρ_y in paper
//! let (ciphertext_y, nonce_y) = key1.encrypt_with_random(
//!     &mut rng,
//!     &(plaintext_y.signed_modulo(key1.n())),
//! )?;
//! // nonce is ρ in paper
//! let (ciphertext_y_by_key1, nonce) = key0.encrypt_with_random(
//!     &mut rng,
//!     &(plaintext_y.signed_modulo(key0.n()))
//! )?;
//! // D in paper
//! let ciphertext_d = key0
//!     .oadd(
//!         &key0.omul(&plaintext_x, &ciphertext_c)?,
//!         &ciphertext_y_by_key1,
//!     )?;
//!
//! // 4. Prover computes a non-interactive proof that plaintext_x and
//! //    plaintext_y are at most `l_x` and `l_y` bits
//!
//! let data = p::Data {
//!     key0,
//!     key1,
//!     c: ciphertext_c,
//!     d: ciphertext_d,
//!     x: ciphertext_x,
//!     y: ciphertext_y,
//! };
//! let pdata = p::PrivateData {
//!     x: plaintext_x,
//!     y: plaintext_y,
//!     nonce,
//!     nonce_y,
//! };
//! let (commitment, proof) =
//!     p::non_interactive::prove(
//!         shared_state_prover,
//!         &aux,
//!         &data,
//!         &pdata,
//!         &security,
//!         &mut rng,
//!     )?;
//!
//! // 5. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<E: Curve>(_: &p::Data<E>, _: &p::Commitment<E>, _: &p::Proof) {  }
//! send(&data, &commitment, &proof);
//!
//! // 6. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//! let r = p::non_interactive::verify(
//!     shared_state_verifier,
//!     &aux,
//!     &data,
//!     &commitment,
//!     &security,
//!     &proof,
//! )?;
//! #
//! # Ok(()) }
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use fast_paillier::{Ciphertext, EncryptionKey, Nonce};
use generic_ec::{Curve, Point};
use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::{Aux, InvalidProof};

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, bit size of +-x
    pub l_x: usize,
    /// l' in paper, bit size of +-y
    pub l_y: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Data<C: Curve> {
    /// N0 in paper, public key that C was encrypted on
    pub key0: EncryptionKey,
    /// N1 in paper, public key that y -> Y was encrypted on
    pub key1: EncryptionKey,
    /// C or C0 in paper, some data encrypted on N0
    pub c: Ciphertext,
    /// D or C in paper, result of affine transformation of C0 with x and y
    pub d: Integer,
    /// Y in paper, y encrypted on N1
    pub y: Ciphertext,
    /// X in paper, obtained as g^x
    pub x: Point<C>,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    /// x or epsilon in paper, preimage of X
    pub x: Integer,
    /// y or delta in paper, preimage of Y
    pub y: Integer,
    /// rho in paper, nonce in encryption of y for additive action
    pub nonce: Nonce,
    /// rho_y in paper, nonce in encryption of y to obtain Y
    pub nonce_y: Nonce,
}

// As described in cggmp21 at page 35
/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Commitment<C: Curve> {
    pub a: Integer,
    pub b_x: Point<C>,
    pub b_y: Integer,
    pub e: Integer,
    pub s: Integer,
    pub f: Integer,
    pub t: Integer,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub beta: Integer,
    pub r: Integer,
    pub r_y: Integer,
    pub gamma: Integer,
    pub m: Integer,
    pub delta: Integer,
    pub mu: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: Integer,
    pub z2: Integer,
    pub z3: Integer,
    pub z4: Integer,
    pub w: Integer,
    pub w_y: Integer,
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use generic_ec::{Curve, Point};
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::common::{fail_if, fail_if_ne, IntegerExt, InvalidProof, InvalidProofReason};
    use crate::Error;

    use super::*;

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        aux: &Aux,
        data: &Data<C>,
        pdata: &PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment<C>, PrivateCommitment), Error> {
        let two_to_l = (Integer::ONE << security.l_x).complete();
        let two_to_l_e = (Integer::ONE << (security.l_x + security.epsilon)).complete();
        let two_to_l_prime_e = (Integer::ONE << (security.l_y + security.epsilon)).complete();
        let hat_n_at_two_to_l_e = (&aux.rsa_modulo * &two_to_l_e).complete();
        let hat_n_at_two_to_l = (&aux.rsa_modulo * &two_to_l).complete();

        let alpha = Integer::from_rng_pm(&two_to_l_e, &mut rng);
        let beta = Integer::from_rng_pm(&two_to_l_prime_e, &mut rng);
        let r = Integer::gen_inversible(data.key0.n(), &mut rng);
        let r_y = Integer::gen_inversible(data.key1.n(), &mut rng);
        let gamma = Integer::from_rng_pm(&hat_n_at_two_to_l_e, &mut rng);
        let delta = Integer::from_rng_pm(&hat_n_at_two_to_l_e, &mut rng);
        let m = Integer::from_rng_pm(&hat_n_at_two_to_l, &mut rng);
        let mu = Integer::from_rng_pm(&hat_n_at_two_to_l, &mut rng);

        let beta_enc_key0 = data.key0.encrypt_with(&beta, &r)?;
        let alpha_at_c = data.key0.omul(&alpha, &data.c)?;
        let a = data.key0.oadd(&alpha_at_c, &beta_enc_key0)?;

        let commitment = Commitment {
            a,
            b_x: Point::<C>::generator() * alpha.to_scalar(),
            b_y: data.key1.encrypt_with(&beta, &r_y)?,
            e: aux.rsa_modulo.combine(&aux.s, &alpha, &aux.t, &gamma)?,
            s: aux.rsa_modulo.combine(&aux.s, &pdata.x, &aux.t, &m)?,
            f: aux.rsa_modulo.combine(&aux.s, &beta, &aux.t, &delta)?,
            t: aux.rsa_modulo.combine(&aux.s, &pdata.y, &aux.t, &mu)?,
        };
        let private_commitment = PrivateCommitment {
            alpha,
            beta,
            r,
            r_y,
            gamma,
            m,
            delta,
            mu,
        };
        Ok((commitment, private_commitment))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<C: Curve>(
        data: &Data<C>,
        pdata: &PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        Ok(Proof {
            z1: (&pcomm.alpha + challenge * &pdata.x).complete(),
            z2: (&pcomm.beta + challenge * &pdata.y).complete(),
            z3: (&pcomm.gamma + challenge * &pcomm.m).complete(),
            z4: (&pcomm.delta + challenge * &pcomm.mu).complete(),
            w: data
                .key0
                .n()
                .combine(&pcomm.r, Integer::ONE, &pdata.nonce, challenge)?,
            w_y: data
                .key1
                .n()
                .combine(&pcomm.r_y, Integer::ONE, &pdata.nonce_y, challenge)?,
        })
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
        // Five equality checks and two range checks
        {
            let lhs = {
                let z1_at_c = data
                    .key0
                    .omul(&proof.z1, &data.c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                let enc = data
                    .key0
                    .encrypt_with(&proof.z2, &proof.w)
                    .map_err(|_| InvalidProofReason::PaillierEnc)?;
                data.key0
                    .oadd(&z1_at_c, &enc)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            let rhs = {
                let e_at_d = data
                    .key0
                    .omul(challenge, &data.d)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key0
                    .oadd(&commitment.a, &e_at_d)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        {
            let lhs = Point::<C>::generator() * proof.z1.to_scalar();
            let rhs = commitment.b_x + data.x * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        {
            let lhs = data
                .key1
                .encrypt_with(&proof.z2, &proof.w_y)
                .map_err(|_| InvalidProofReason::PaillierEnc)?;
            let rhs = {
                let e_at_y = data
                    .key1
                    .omul(challenge, &data.y)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key1
                    .oadd(&commitment.b_y, &e_at_y)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(3), lhs, rhs)?;
        }
        fail_if_ne(
            InvalidProofReason::EqualityCheck(4),
            aux.rsa_modulo
                .combine(&aux.s, &proof.z1, &aux.t, &proof.z3)?,
            aux.rsa_modulo
                .combine(&commitment.e, Integer::ONE, &commitment.s, challenge)?,
        )?;
        fail_if_ne(
            InvalidProofReason::EqualityCheck(5),
            aux.rsa_modulo
                .combine(&aux.s, &proof.z2, &aux.t, &proof.z4)?,
            aux.rsa_modulo
                .combine(&commitment.f, Integer::ONE, &commitment.t, challenge)?,
        )?;
        fail_if(
            InvalidProofReason::RangeCheck(6),
            proof
                .z1
                .is_in_pm(&(Integer::ONE << (security.l_x + security.epsilon)).complete()),
        )?;
        fail_if(
            InvalidProofReason::RangeCheck(7),
            proof
                .z2
                .is_in_pm(&(Integer::ONE << (security.l_y + security.epsilon)).complete()),
        )?;
        Ok(())
    }

    /// Generate random challenge
    pub fn challenge<R>(security: &SecurityParams, rng: &mut R) -> Integer
    where
        R: RngCore,
    {
        Integer::from_rng_pm(&security.q, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use digest::{typenum::U32, Digest};
    use generic_ec::{hash_to_curve::FromHash, Curve, Scalar};
    use rand_core::RngCore;

    use crate::{Error, InvalidProof};

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<C: Curve, R: RngCore, D>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        pdata: &PrivateData,
        security: &SecurityParams,
        rng: R,
    ) -> Result<(Commitment<C>, Proof), Error>
    where
        Scalar<C>: FromHash,
        D: Digest<OutputSize = U32>,
    {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge)?;
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

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        commitment: &Commitment<C>,
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
                .chain_update((security.l_x as u64).to_le_bytes())
                .chain_update((security.l_y as u64).to_le_bytes())
                .chain_update((security.epsilon as u64).to_le_bytes())
                .chain_update(data.key0.n().to_digits::<u8>(order))
                .chain_update(data.key1.n().to_digits::<u8>(order))
                .chain_update(data.c.to_digits::<u8>(order))
                .chain_update(data.d.to_digits::<u8>(order))
                .chain_update(data.y.to_digits::<u8>(order))
                .chain_update(data.x.to_bytes(true))
                .chain_update(commitment.a.to_digits::<u8>(order))
                .chain_update(commitment.b_x.to_bytes(true))
                .chain_update(commitment.b_y.to_digits::<u8>(order))
                .chain_update(commitment.e.to_digits::<u8>(order))
                .chain_update(commitment.s.to_digits::<u8>(order))
                .chain_update(commitment.f.to_digits::<u8>(order))
                .chain_update(commitment.t.to_digits::<u8>(order))
                .finalize()
        };
        let mut rng = crate::common::rng::HashRng::new(hash);
        super::interactive::challenge(security, &mut rng)
    }
}

#[cfg(test)]
mod test {
    use generic_ec::Point;
    use generic_ec::{hash_to_curve::FromHash, Curve, Scalar};
    use rug::{Complete, Integer};

    use crate::common::test::random_key;
    use crate::common::{IntegerExt, InvalidProofReason};

    fn run<R: rand_core::RngCore + rand_core::CryptoRng, C: Curve>(
        rng: &mut R,
        security: super::SecurityParams,
        x: Integer,
        y: Integer,
    ) -> Result<(), crate::common::InvalidProof>
    where
        Scalar<C>: FromHash,
    {
        let dk0 = random_key(rng).unwrap();
        let dk1 = random_key(rng).unwrap();
        let ek0 = dk0.encryption_key();
        let ek1 = dk1.encryption_key();

        let (c, _) = {
            let plaintext = Integer::from_rng_pm(ek0.half_n(), rng);
            ek0.encrypt_with_random(rng, &plaintext).unwrap()
        };

        let (y_enc_ek1, rho_y) = ek1.encrypt_with_random(rng, &y).unwrap();

        let (y_enc_ek0, rho) = ek0.encrypt_with_random(rng, &y).unwrap();
        let x_at_c = ek0.omul(&x, &c).unwrap();
        let d = ek0.oadd(&x_at_c, &y_enc_ek0).unwrap();

        let data = super::Data {
            key0: ek0,
            key1: ek1,
            c,
            d,
            y: y_enc_ek1,
            x: x.to_scalar::<C>() * Point::generator(),
        };
        let pdata = super::PrivateData {
            x,
            y,
            nonce: rho,
            nonce_y: rho_y,
        };

        let aux = crate::common::test::aux(rng);

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
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l_x: 1024,
            l_y: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).into(),
        };
        let x = Integer::from_rng_pm(&(Integer::ONE << security.l_x).complete(), &mut rng);
        let y = Integer::from_rng_pm(&(Integer::ONE << security.l_y).complete(), &mut rng);
        run::<_, C>(&mut rng, security, x, y).expect("proof failed");
    }

    fn failing_on_additive<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l_x: 1024,
            l_y: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let x = Integer::from_rng_pm(&(Integer::ONE << security.l_x).complete(), &mut rng);
        let y = (Integer::ONE << (security.l_y + security.epsilon)).complete() + 1;
        let r = run::<_, C>(&mut rng, security, x, y).expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::RangeCheck(7) => (),
            e => panic!("proof should not fail with: {e:?}"),
        }
    }

    fn failing_on_multiplicative<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l_x: 1024,
            l_y: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let x = (Integer::ONE << (security.l_x + security.epsilon)).complete() + 1;
        let y = Integer::from_rng_pm(&(Integer::ONE << security.l_y).complete(), &mut rng);
        let r = run::<_, C>(&mut rng, security, x, y).expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::RangeCheck(6) => (),
            e => panic!("proof should not fail with: {e:?}"),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec::curves::Secp256r1>()
    }
    #[test]
    fn failing_p256_add() {
        failing_on_additive::<generic_ec::curves::Secp256r1>()
    }
    #[test]
    fn failing_p256_mul() {
        failing_on_multiplicative::<generic_ec::curves::Secp256r1>()
    }

    #[test]
    fn passing_million() {
        passing_test::<crate::curve::C>()
    }
    #[test]
    fn failing_million_add() {
        failing_on_additive::<crate::curve::C>()
    }
    #[test]
    fn failing_million_mul() {
        failing_on_multiplicative::<crate::curve::C>()
    }
}
