//! ZK-proof of knowledge of plaintext of paillier encryption. Called Пdec or
//! Rdec in the CGGMP21 paper.
//!
//! ## Description
//! A party P has `key`, `pkey` - public and private keys in paillier
//! cryptosystem. P has plaintext y, `nonce` and `c = key.encrypt(y, nonce)`.
//! P also prepares a number q and x such that `x = y mod q`.
//! P wants to prove that it knows y while disclosing x, c and key.
//!
//! ## Example
//!
//! ```no_run
//! # use paillier_zk::unknown_order::BigNumber;
//! use paillier_zk::paillier_decryption_modulo_q as p;
//! use generic_ec::hash_to_curve::Tag;
//! let shared_state_prover = sha2::Sha256::default();
//! let shared_state_verifier = sha2::Sha256::default();
//! let mut rng = rand_core::OsRng::default();
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
//! let security_prime = BigNumber::prime(256);
//!
//! let aux = p::Aux { s, t, rsa_modulo };
//!
//! // 1. Setup: prover prepares the paillier keys and q
//!
//! let private_key = libpaillier::DecryptionKey::random().unwrap();
//! let key = libpaillier::EncryptionKey::from(&private_key);
//! let q = BigNumber::prime(1024);
//!
//! // 2. Setup: prover has some plaintext and encrypts it
//!
//! let plaintext: BigNumber = 228.into();
//! let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
//! let random_koef = BigNumber::from_rng(&q, &mut rng);
//! let x = &plaintext + &q * random_koef;
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 128,
//! };
//!
//! let data = p::Data { key, c: ciphertext, q, x };
//! let pdata = p::PrivateData { y: plaintext, nonce };
//! let (commitment, proof) =
//!     p::non_interactive::prove(shared_state_prover, &aux, &data, &pdata, &security, &mut rng)
//!         .expect("could not compute proof");
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Proof) { todo!() }
//! # fn recv() -> (p::Data, p::Commitment, p::Proof) { todo!() }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! let (data, commitment, proof) = recv();
//! p::non_interactive::verify(shared_state_verifier, &aux, &data, &commitment, &proof);
//!
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover
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
    /// l in paper, bit size of plaintext
    pub l: usize,
    /// Epsilon in paper, range extension and security parameter for x
    pub epsilon: usize,
}

/// Public data that both parties know
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Data {
    /// q in paper, modulo value such that `x = y mod q`
    pub q: BigNumber,
    /// N0 in paper, public key that y -> C was encrypted on
    pub key: EncryptionKey,
    /// C in paper
    pub c: Ciphertext,
    /// x in paper
    pub x: BigNumber,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    /// y in paper, plaintext value of C
    pub y: BigNumber,
    /// rho in paper, nonce of encryption y -> C
    pub nonce: Nonce,
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// Prover's first message, obtained by [`interactive::commit`]
pub struct Commitment {
    pub s: BigNumber,
    pub t: BigNumber,
    pub a: Ciphertext,
    pub gamma: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub mu: BigNumber,
    pub nu: BigNumber,
    pub r: Nonce,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = BigNumber;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]. Consists of M proofs for each challenge
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub w: BigNumber,
}

pub use crate::common::Aux;

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use crate::unknown_order::BigNumber;
    use rand_core::RngCore;

    use crate::common::{combine, InvalidProof, ProtocolError};

    use super::{
        Aux, Challenge, Commitment, Data, PrivateCommitment, PrivateData, Proof, SecurityParams,
    };

    /// Create random commitment
    pub fn commit<R: RngCore>(
        aux: &Aux,
        data: &Data,
        pdata: &PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment, PrivateCommitment), ProtocolError> {
        let two_to_l_e = BigNumber::one() << (security.l + security.epsilon + 1);
        let modulo_l = (BigNumber::one() << (security.l + 1)) * &aux.rsa_modulo;
        let modulo_l_e = &two_to_l_e * &aux.rsa_modulo;

        let alpha = BigNumber::from_rng(&two_to_l_e, &mut rng);
        let mu = BigNumber::from_rng(&modulo_l, &mut rng);
        let nu = BigNumber::from_rng(&modulo_l_e, &mut rng);

        let (a, r) = data
            .key
            .encrypt(alpha.to_bytes(), None)
            .ok_or(ProtocolError::EncryptionFailed)?;

        let commitment = Commitment {
            s: combine(&aux.s, &pdata.y, &aux.t, &mu, &aux.rsa_modulo),
            t: combine(&aux.s, &alpha, &aux.t, &nu, &aux.rsa_modulo),
            a,
            gamma: &alpha % &data.q,
        };
        let private_commitment = PrivateCommitment { alpha, mu, nu, r };
        Ok((commitment, private_commitment))
    }

    /// Generate random challenge
    pub fn challenge<R: RngCore>(data: &Data, rng: &mut R) -> Challenge {
        // double the range to account for +-
        let m = BigNumber::from(2) * &data.q;
        BigNumber::from_rng(&m, rng)
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        data: &Data,
        pdata: &PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Proof {
        Proof {
            z1: &pcomm.alpha + challenge * &pdata.y,
            z2: &pcomm.nu + challenge * &pcomm.mu,
            w: combine(
                &pcomm.r,
                &BigNumber::one(),
                &pdata.nonce,
                challenge,
                data.key.n(),
            ),
        }
    }

    /// Verify the proof
    pub fn verify(
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
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
        // Three equality checks
        {
            let (lhs, _) = data
                .key
                .encrypt(proof.z1.to_bytes(), Some(proof.w.clone()))
                .ok_or(InvalidProof::EncryptionFailed)?;
            let rhs = combine(&commitment.a, &one, &data.c, challenge, data.key.nn());
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(1))?;
        }
        {
            let lhs = &proof.z1 % &data.q;
            let rhs = commitment
                .gamma
                .modadd(&challenge.modmul(&data.x, &data.q), &data.q);
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(2))?;
        }
        {
            let lhs = combine(&aux.s, &proof.z1, &aux.t, &proof.z2, &aux.rsa_modulo);
            let rhs = combine(
                &commitment.t,
                &one,
                &commitment.s,
                challenge,
                &aux.rsa_modulo,
            );
            fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(3))?;
        }

        Ok(())
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use crate::unknown_order::BigNumber;
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use crate::common::{InvalidProof, ProtocolError};

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<R: RngCore, D>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        pdata: &PrivateData,
        security: &SecurityParams,
        rng: R,
    ) -> Result<(Commitment, Proof), ProtocolError>
    where
        D: Digest<OutputSize = U32>,
    {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge(shared_state, aux, data, &comm);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge);
        Ok((comm, proof))
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<D>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
    ) -> Challenge
    where
        D: Digest<OutputSize = U32>,
    {
        use rand_core::SeedableRng;
        let seed = shared_state
            .chain_update(aux.s.to_bytes())
            .chain_update(aux.t.to_bytes())
            .chain_update(aux.rsa_modulo.to_bytes())
            .chain_update(data.q.to_bytes())
            .chain_update(data.key.to_bytes())
            .chain_update(data.c.to_bytes())
            .chain_update(data.x.to_bytes())
            .chain_update(commitment.s.to_bytes())
            .chain_update(commitment.t.to_bytes())
            .chain_update(commitment.a.to_bytes())
            .chain_update(commitment.gamma.to_bytes())
            .finalize();
        let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed.into());
        let m = BigNumber::from(2) * &data.q;
        BigNumber::from_rng(&m, &mut rng)
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<D>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        proof: &Proof,
    ) -> Result<(), InvalidProof>
    where
        D: Digest<OutputSize = U32>,
    {
        let challenge = challenge(shared_state, aux, data, commitment);
        super::interactive::verify(aux, data, commitment, &challenge, proof)
    }
}

#[cfg(test)]
mod test {
    use libpaillier::unknown_order::BigNumber;

    #[test]
    fn passing_test() {
        let mut rng = rand_core::OsRng::default();

        let aux = crate::common::test::aux(&mut rng);

        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
        };
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(28);
        let hiddentext = BigNumber::from(228);
        let modulo = BigNumber::from(100);
        assert_eq!(&plaintext % &modulo, &hiddentext % &modulo);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();

        let data = super::Data {
            q: modulo,
            key: key0,
            c: ciphertext,
            x: hiddentext,
        };
        let pdata = super::PrivateData {
            y: plaintext,
            nonce,
        };

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
        let r = super::non_interactive::verify(shared_state, &aux, &data, &commitment, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }

    #[test]
    fn failing_wrong_hidden() {
        let mut rng = rand_core::OsRng::default();

        let aux = crate::common::test::aux(&mut rng);

        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
        };
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(28);
        let hiddentext = BigNumber::from(322);
        let modulo = BigNumber::from(100);
        assert_ne!(&plaintext % &modulo, &hiddentext % &modulo);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();

        let data = super::Data {
            q: modulo,
            key: key0,
            c: ciphertext,
            x: hiddentext,
        };
        let pdata = super::PrivateData {
            y: plaintext,
            nonce,
        };

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
        let r = super::non_interactive::verify(shared_state, &aux, &data, &commitment, &proof);
        if r.is_ok() {
            panic!("proof should not pass");
        }
    }

    #[test]
    fn failing_wrong_plain() {
        let mut rng = rand_core::OsRng::default();

        let aux = crate::common::test::aux(&mut rng);

        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
        };
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let wrong_plaintext = BigNumber::from(42);
        let plaintext = BigNumber::from(28);
        let hiddentext = BigNumber::from(228);
        let modulo = BigNumber::from(100);
        let (ciphertext, nonce) = key0.encrypt(wrong_plaintext.to_bytes(), None).unwrap();

        let data = super::Data {
            q: modulo,
            key: key0,
            c: ciphertext,
            x: hiddentext,
        };
        let pdata = super::PrivateData {
            y: plaintext,
            nonce,
        };

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
        let r = super::non_interactive::verify(shared_state, &aux, &data, &commitment, &proof);
        if r.is_ok() {
            panic!("proof should not pass");
        }
    }

    // Following motivation outlined in
    // [crate::paillier_encryption_in_range::test::rejected_with_probability_1_over_2],
    // I would like to make a similar borderline test, but no security estimate
    // was given in the paper and this proof differs significantly from others
    // in this library, so I have to omit the test.
}
