//! ZK-proof of paillier operation with group commitment in range. Called Пaff-g
//! or Raff-g in the CGGMP21 paper.
//!
//! ## Description
//!
//! A party P performs a paillier affine operation with C, Y, and X
//! obtaining `D = C*X + Y`. `X` and `Y` are encrypted values of `x` and `y`. P
//! then wants to prove that `y` and `x` are at most `L+1` and `L'+1` bits,
//! correspondingly, and P doesn't want to disclose none of the plaintexts
//!
//! Given:
//! - `key0`, `pkey0`, `key1`, `pkey1` - pairs of public and private keys in
//!   paillier cryptosystem
//! - `nonce_c`, `nonce_y`, `nonce` - nonces in paillier encryption
//! - `c`, `x`, `y` - some numbers
//! - `q`, `g` such that `<g> = Zq*` - prime order group
//! - `C = key0.encrypt(c, nonce_c)`
//! - `Y' = key0.encrypt(y, nonce)`
//! - `Y = key1.encrypt(y, nonce_y)`
//! - `X = g * x`
//! - `D = key0.affine_operation!{ X * C + Y' }`, i.e.
//!   `pkey0.decrypt(D) = x * c + y`
//!
//! Prove:
//! - `bitsize(y) <= L`
//! - `bitsize(x) <= L'`
//!
//! Disclosing only: `key0`, `key1`, `C`, `D`, `Y`, `X`
//!
//! ## Example
//!
//! ``` no_run
//! # use paillier_zk::unknown_order::BigNumber;
//! use paillier_zk::paillier_affine_operation_in_range as p;
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
//! // this key is used to decrypt C and D, and also Y in the affine operation
//! let private_key0 = libpaillier::DecryptionKey::random().unwrap();
//! let key0 = libpaillier::EncryptionKey::from(&private_key0);
//! // this key is used to decrypt Y in this ZK-protocol
//! let private_key1 = libpaillier::DecryptionKey::random().unwrap();
//! let key1 = libpaillier::EncryptionKey::from(&private_key1);
//!
//! // 2. Setup: prover prepares the group used to encrypt x
//!
//! type C = generic_ec_curves::Secp256r1;
//! let g = generic_ec::Point::<C>::generator();
//!
//! // 3. Setup: prover prepares all plain texts
//!
//! // c in paper
//! let plaintext_orig = BigNumber::from(100);
//! // x in paper
//! let plaintext_mult = BigNumber::from(2);
//! // y in paper
//! let plaintext_add = BigNumber::from(28);
//!
//! // 4. Setup: prover encrypts everything on correct keys and remembers some nonces
//!
//! // C in paper
//! let (ciphertext_orig, _) = key0.encrypt(plaintext_orig.to_bytes(), None).unwrap();
//! // X in paper
//! let ciphertext_mult = g * paillier_zk::convert_scalar(&plaintext_mult);
//! // Y' in further docs, and ρy in paper
//! let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), None).unwrap();
//! // Y and ρ in paper
//! let (ciphertext_add_action, nonce) = key0.encrypt(plaintext_add.to_bytes(), None).unwrap();
//! // D in paper
//! let transformed = key0
//!     .add(
//!         &key0.mul(&ciphertext_orig, &plaintext_mult).unwrap(),
//!         &ciphertext_add_action,
//!     )
//!     .unwrap();
//!
//! // 5. Prover computes a non-interactive proof that plaintext_add and
//! //    plaintext_mult are at most L and L' bits
//!
//! let mut rng = rand_core::OsRng::default();
//!
//! let security = p::SecurityParams {
//!     l_x: 1024,
//!     l_y: 1024,
//!     epsilon: 128,
//!     q: BigNumber::prime_from_rng(128, &mut rng),
//! };
//!
//! let data = p::Data {
//!     key0,
//!     key1,
//!     c: ciphertext_orig,
//!     d: transformed,
//!     y: ciphertext_add,
//!     x: ciphertext_mult,
//! };
//! let pdata = p::PrivateData {
//!     x: plaintext_mult,
//!     y: plaintext_add,
//!     nonce,
//!     nonce_y,
//! };
//! let (commitment, proof) =
//!     p::non_interactive::prove(shared_state_prover, &aux, &data, &pdata, &security, rng).expect("proof failed");
//!
//! // 6. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<C: Curve>(_: &p::Data<C>, _: &p::Commitment<C>, _: &p::Proof) { todo!() }
//! # fn recv<C: Curve>() -> (p::Data<C>, p::Commitment<C>, p::Proof) { todo!() }
//! send(&data, &commitment, &proof);
//!
//! // 7. Verifier receives the data and the proof and verifies it
//!
//! let (data, commitment, proof) = recv::<C>();
//! let r = p::non_interactive::verify(shared_state_verifier, &aux, &data, &commitment, &security, &proof);
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use generic_ec::{Curve, Point};
use libpaillier::{Ciphertext, EncryptionKey, Nonce};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::InvalidProof;
use crate::unknown_order::BigNumber;

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
    pub q: BigNumber,
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
    pub d: BigNumber,
    /// Y in paper, y encrypted on N1
    pub y: Ciphertext,
    /// X in paper, obtained as g^x
    pub x: Point<C>,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    /// x or epsilon in paper, preimage of X
    pub x: BigNumber,
    /// y or delta in paper, preimage of Y
    pub y: BigNumber,
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
    pub a: BigNumber,
    pub b_x: Point<C>,
    pub b_y: BigNumber,
    pub e: BigNumber,
    pub s: BigNumber,
    pub f: BigNumber,
    pub t: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub beta: BigNumber,
    pub r: BigNumber,
    pub r_y: BigNumber,
    pub gamma: BigNumber,
    pub m: BigNumber,
    pub delta: BigNumber,
    pub mu: BigNumber,
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
    pub z4: BigNumber,
    pub w: BigNumber,
    pub w_y: BigNumber,
}

pub use crate::common::Aux;

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {

    use crate::unknown_order::BigNumber;
    use generic_ec::{Curve, Point};
    use rand_core::RngCore;

    use crate::common::{combine, convert_scalar, gen_inversible, InvalidProof, ProtocolError};

    use super::*;

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        aux: &Aux,
        data: &Data<C>,
        pdata: &PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment<C>, PrivateCommitment), ProtocolError> {
        let two_to_l = BigNumber::one() << (security.l_x + 1);
        let two_to_l_e = BigNumber::one() << (security.l_x + security.epsilon + 1);
        let modulo_l = two_to_l * &aux.rsa_modulo;
        let modulo_l_e = &two_to_l_e * &aux.rsa_modulo;

        let alpha = BigNumber::from_rng(&two_to_l_e, &mut rng);
        let beta = BigNumber::from_rng(
            &(BigNumber::one() << (security.l_y + security.epsilon + 1)),
            &mut rng,
        );
        let r = gen_inversible(data.key0.n(), &mut rng);
        let r_y = gen_inversible(data.key1.n(), &mut rng);
        let gamma = BigNumber::from_rng(&modulo_l_e, &mut rng);
        let m = BigNumber::from_rng(&modulo_l, &mut rng);
        let delta = BigNumber::from_rng(&modulo_l_e, &mut rng);
        let mu = BigNumber::from_rng(&modulo_l, &mut rng);

        let a_add = data
            .key0
            .encrypt(beta.to_bytes(), Some(r.clone()))
            .ok_or(ProtocolError::EncryptionFailed)?
            .0;
        let c_to_alpha = data
            .key0
            .mul(&data.c, &alpha)
            .ok_or(ProtocolError::EncryptionFailed)?;
        let a = data
            .key0
            .add(&c_to_alpha, &a_add)
            .ok_or(ProtocolError::EncryptionFailed)?;
        let commitment = Commitment {
            a,
            b_x: Point::<C>::generator() * convert_scalar(&alpha),
            b_y: data
                .key1
                .encrypt(beta.to_bytes(), Some(r_y.clone()))
                .ok_or(ProtocolError::EncryptionFailed)?
                .0,
            e: combine(&aux.s, &alpha, &aux.t, &gamma, &aux.rsa_modulo),
            s: combine(&aux.s, &pdata.x, &aux.t, &m, &aux.rsa_modulo),
            f: combine(&aux.s, &beta, &aux.t, &delta, &aux.rsa_modulo),
            t: combine(&aux.s, &pdata.y, &aux.t, &mu, &aux.rsa_modulo),
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
    ) -> Proof {
        Proof {
            z1: &pcomm.alpha + challenge * &pdata.x,
            z2: &pcomm.beta + challenge * &pdata.y,
            z3: &pcomm.gamma + challenge * &pcomm.m,
            z4: &pcomm.delta + challenge * &pcomm.mu,
            w: combine(
                &pcomm.r,
                &BigNumber::one(),
                &pdata.nonce,
                challenge,
                data.key0.n(),
            ),
            w_y: combine(
                &pcomm.r_y,
                &BigNumber::one(),
                &pdata.nonce_y,
                challenge,
                data.key1.n(),
            ),
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
        fn fail_if(msg: InvalidProof, b: bool) -> Result<(), InvalidProof> {
            if b {
                Ok(())
            } else {
                Err(msg)
            }
        }
        // Five equality checks and two range checks
        {
            let enc = data
                .key0
                .encrypt(proof.z2.to_bytes(), Some(proof.w.clone()))
                .ok_or(InvalidProof::EncryptionFailed)?
                .0;
            let lhs = data
                .key0
                .add(
                    &data
                        .key0
                        .mul(&data.c, &proof.z1)
                        .ok_or(InvalidProof::EncryptionFailed)?,
                    &enc,
                )
                .ok_or(InvalidProof::EncryptionFailed)?;
            let rhs = combine(&commitment.a, &one, &data.d, challenge, data.key0.nn());
            fail_if(InvalidProof::EqualityCheckFailed(1), lhs == rhs)?;
        }
        {
            let lhs = Point::<C>::generator() * convert_scalar(&proof.z1);
            let rhs = commitment.b_x + data.x * convert_scalar(challenge);
            fail_if(InvalidProof::EqualityCheckFailed(2), lhs == rhs)?;
        }
        {
            let lhs = data
                .key1
                .encrypt(proof.z2.to_bytes(), Some(proof.w_y.clone()))
                .ok_or(InvalidProof::EncryptionFailed)?
                .0;
            let rhs = combine(&commitment.b_y, &one, &data.y, challenge, data.key1.nn());
            fail_if(InvalidProof::EqualityCheckFailed(3), lhs == rhs)?;
        }
        fail_if(
            InvalidProof::EqualityCheckFailed(4),
            combine(&aux.s, &proof.z1, &aux.t, &proof.z3, &aux.rsa_modulo)
                == combine(
                    &commitment.e,
                    &one,
                    &commitment.s,
                    challenge,
                    &aux.rsa_modulo,
                ),
        )?;
        fail_if(
            InvalidProof::EqualityCheckFailed(5),
            combine(&aux.s, &proof.z2, &aux.t, &proof.z4, &aux.rsa_modulo)
                == combine(
                    &commitment.f,
                    &one,
                    &commitment.t,
                    challenge,
                    &aux.rsa_modulo,
                ),
        )?;
        fail_if(
            InvalidProof::RangeCheckFailed(6),
            proof.z1 <= &one << (security.l_x + security.epsilon + 1),
        )?;
        fail_if(
            InvalidProof::RangeCheckFailed(7),
            proof.z2 <= &one << (security.l_y + security.epsilon + 1),
        )?;
        Ok(())
    }

    /// Generate random challenge
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

    use crate::unknown_order::BigNumber;
    use generic_ec::{hash_to_curve::FromHash, Curve, Scalar};
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use crate::common::{InvalidProof, ProtocolError};

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
    ) -> Result<(Commitment<C>, Proof), ProtocolError>
    where
        Scalar<C>: FromHash,
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

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: D,
        aux: &Aux,
        data: &Data<C>,
        commitment: &Commitment<C>,
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
            .chain_update(data.key0.to_bytes())
            .chain_update(data.key1.to_bytes())
            .chain_update(data.c.to_bytes())
            .chain_update(data.d.to_bytes())
            .chain_update(data.y.to_bytes())
            .chain_update(data.x.to_bytes(true))
            .chain_update(commitment.a.to_bytes())
            .chain_update(commitment.b_x.to_bytes(true))
            .chain_update(commitment.b_y.to_bytes())
            .chain_update(commitment.e.to_bytes())
            .chain_update(commitment.s.to_bytes())
            .chain_update(commitment.f.to_bytes())
            .chain_update(commitment.t.to_bytes())
            .chain_update((security.l_x as u64).to_le_bytes())
            .chain_update((security.l_y as u64).to_le_bytes())
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

    use crate::common::convert_scalar;
    use crate::unknown_order::BigNumber;

    fn random_key<R: rand_core::RngCore>(rng: &mut R) -> Option<libpaillier::DecryptionKey> {
        let p = BigNumber::prime_from_rng(1024, rng);
        let q = BigNumber::prime_from_rng(1024, rng);
        libpaillier::DecryptionKey::with_primes_unchecked(&p, &q)
    }

    fn nonce<R: rand_core::RngCore>(rng: &mut R, n: &BigNumber) -> Option<BigNumber> {
        Some(BigNumber::from_rng(n, rng))
    }

    fn run<R: rand_core::RngCore, C: Curve>(
        mut rng: R,
        security: super::SecurityParams,
        plaintext_orig: BigNumber,
        plaintext_mult: BigNumber,
        plaintext_add: BigNumber,
    ) -> Result<(), crate::common::InvalidProof>
    where
        Scalar<C>: FromHash,
    {
        let affined = &plaintext_mult * &plaintext_orig + &plaintext_add;

        let private_key0 = random_key(&mut rng).unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);
        let private_key1 = random_key(&mut rng).unwrap();
        let key1 = libpaillier::EncryptionKey::from(&private_key1);
        let g = generic_ec::Point::<C>::generator();
        let (ciphertext, _) = key0
            .encrypt(affined.to_bytes(), nonce(&mut rng, key0.n()))
            .unwrap();
        let (ciphertext_orig, _) = key0
            .encrypt(plaintext_orig.to_bytes(), nonce(&mut rng, key0.n()))
            .unwrap();
        let ciphertext_mult = g * convert_scalar(&plaintext_mult);
        let nonce_y = nonce(&mut rng, key1.n());
        let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), nonce_y).unwrap();
        let nonce = nonce(&mut rng, key0.n());
        let (ciphertext_add_action, nonce) = key0.encrypt(plaintext_add.to_bytes(), nonce).unwrap();
        // verify that D is obtained from affine transformation of C
        let transformed = key0
            .add(
                &key0.mul(&ciphertext_orig, &plaintext_mult).unwrap(),
                &ciphertext_add_action,
            )
            .unwrap();
        assert_eq!(
            private_key0.decrypt(&transformed).unwrap(),
            private_key0.decrypt(&ciphertext).unwrap(),
        );
        let data = super::Data {
            key0,
            key1,
            c: ciphertext_orig,
            d: transformed,
            y: ciphertext_add,
            x: ciphertext_mult,
        };
        let pdata = super::PrivateData {
            x: plaintext_mult,
            y: plaintext_add,
            nonce,
            nonce_y,
        };

        let p = BigNumber::prime_from_rng(1024, &mut rng);
        let q = BigNumber::prime_from_rng(1024, &mut rng);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

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
            l_x: 1024,
            l_y: 1024,
            epsilon: 256,
            q: BigNumber::prime_from_rng(128, &mut rng),
        };
        let plaintext_orig = BigNumber::from(100);
        let plaintext_mult = BigNumber::from(1) << (security.l_x + 1);
        let plaintext_add = BigNumber::from(1) << (security.l_y + 1);
        let r = run(rng, security, plaintext_orig, plaintext_mult, plaintext_add);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }

    fn failing_on_additive<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_core::OsRng::default();
        let security = super::SecurityParams {
            l_x: 1024,
            l_y: 1024,
            epsilon: 256,
            q: BigNumber::prime_from_rng(128, &mut rng),
        };
        let plaintext_orig = BigNumber::from(100);
        let plaintext_mult = BigNumber::from(1) << (security.l_x + 1);
        let plaintext_add = (BigNumber::from(1) << (security.l_y + security.epsilon + 1)) + 1;
        let r = run(rng, security, plaintext_orig, plaintext_mult, plaintext_add);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(crate::common::InvalidProof::RangeCheckFailed(7)) => (),
            Err(e) => panic!("proof should not fail with: {:?}", e),
        }
    }

    fn failing_on_multiplicative<C: Curve>()
    where
        Scalar<C>: FromHash,
    {
        let mut rng = rand_core::OsRng::default();
        let security = super::SecurityParams {
            l_x: 1024,
            l_y: 1024,
            epsilon: 256,
            q: BigNumber::prime_from_rng(128, &mut rng),
        };
        let plaintext_orig = BigNumber::from(100);
        let plaintext_mult = (BigNumber::from(1) << (security.l_x + security.epsilon + 1)) + 1;
        let plaintext_add = BigNumber::from(1) << (security.l_y + 1);
        let r = run(rng, security, plaintext_orig, plaintext_mult, plaintext_add);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(crate::common::InvalidProof::RangeCheckFailed(6)) => (),
            Err(e) => panic!("proof should not fail with: {:?}", e),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec_curves::rust_crypto::Secp256r1>()
    }
    #[test]
    fn failing_p256_add() {
        failing_on_additive::<generic_ec_curves::rust_crypto::Secp256r1>()
    }
    #[test]
    fn failing_p256_mul() {
        failing_on_multiplicative::<generic_ec_curves::rust_crypto::Secp256r1>()
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

    // see notes in
    // [crate::paillier_encryption_in_range::test::rejected_with_probability_1_over_2]
    // for motivation and design of the following two tests

    #[test]
    fn mul_rejected_with_probability_1_over_2() {
        use rand_core::SeedableRng;
        fn maybe_rejected(mut rng: rand_chacha::ChaCha20Rng) -> bool {
            let security = super::SecurityParams {
                l_x: 1024,
                l_y: 1024,
                epsilon: 130,
                q: BigNumber::prime_from_rng(128, &mut rng),
            };
            let plaintext_orig = BigNumber::from(100);
            let plaintext_mult = (BigNumber::from(1) << (security.l_x + 1)) - 1;
            let plaintext_add = BigNumber::from(1) << (security.l_y / 2);
            let r = run::<_, generic_ec_curves::rust_crypto::Secp256r1>(
                rng,
                security,
                plaintext_orig,
                plaintext_mult,
                plaintext_add,
            );
            match r {
                Ok(()) => true,
                Err(crate::common::InvalidProof::RangeCheckFailed(6)) => false,
                Err(e) => panic!("proof should not fail with: {:?}", e),
            }
        }

        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(0);
        assert!(!maybe_rejected(rng), "should fail");
        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(1);
        assert!(maybe_rejected(rng), "should pass");
    }

    #[test]
    fn add_rejected_with_probability_1_over_2() {
        use rand_core::SeedableRng;
        fn maybe_rejected(mut rng: rand_chacha::ChaCha20Rng) -> bool {
            let security = super::SecurityParams {
                l_x: 1024,
                l_y: 1024,
                epsilon: 130,
                q: BigNumber::prime_from_rng(128, &mut rng),
            };
            let plaintext_orig = BigNumber::from(100);
            let plaintext_mult = BigNumber::from(1) << (security.l_x / 2);
            let plaintext_add = (BigNumber::from(1) << (security.l_y + 1)) + 1;
            let r = run::<_, generic_ec_curves::rust_crypto::Secp256r1>(
                rng,
                security,
                plaintext_orig,
                plaintext_mult,
                plaintext_add,
            );
            match r {
                Ok(()) => true,
                Err(crate::common::InvalidProof::RangeCheckFailed(7)) => false,
                Err(e) => panic!("proof should not fail with: {:?}", e),
            }
        }

        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(0);
        assert!(maybe_rejected(rng), "should pass");
        let rng = rand_chacha::ChaCha20Rng::seed_from_u64(1);
        assert!(!maybe_rejected(rng), "should fail");
    }
}
