//! ZK-proof of paillier encryption in range. Called Пenc or Renc in the CGGMP21
//! paper.
//!
//! ## Description
//!
//! A party P has `key`, `pkey` - public and private keys in paillier
//! cryptosystem. P also has `plaintext`, `nonce`, and
//! `ciphertext = key.encrypt(plaintext, nonce)`.
//!
//! P wants to prove that `plaintext` is at most `L` bits, without disclosing
//! it, the `pkey`, and `nonce`

//! ## Example
//!
//! ``` no_run
//! # use paillier_zk::unknown_order::BigNumber;
//! use paillier_zk::paillier_encryption_in_range as p;
//! use generic_ec::hash_to_curve::Tag;
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
//! let security_prime = BigNumber::prime(256);
//!
//! let aux = p::Aux { s, t, rsa_modulo };
//!
//! // 1. Setup: prover prepares the paillier keys
//!
//! let private_key = libpaillier::DecryptionKey::random().unwrap();
//! let key = libpaillier::EncryptionKey::from(&private_key);
//!
//! // 2. Setup: prover has some plaintext and encrypts it
//!
//! let plaintext: BigNumber = 228.into();
//! let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 128,
//!     q: security_prime,
//! };
//!
//! let rng = rand_core::OsRng::default();
//! let data = p::Data { key, ciphertext };
//! let pdata = p::PrivateData { plaintext, nonce };
//! let (commitment, proof) =
//!     p::non_interactive::prove(shared_state_prover, &aux, &data, &pdata, &security, rng);
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
//! p::non_interactive::verify(shared_state_verifier, &aux, &data, &commitment, &security, &proof);
//!
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use crate::unknown_order::BigNumber;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};

pub use crate::common::InvalidProof;

#[derive(Debug, Clone)]
pub struct SecurityParams {
    /// l in paper, bit size of plaintext
    pub l: usize,
    /// Epsilon in paper, range extension and security parameter for x
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: BigNumber,
}

/// Public data that both parties know
#[derive(Debug, Clone)]
pub struct Data {
    /// N0 in paper, public key that k -> K was encrypted on
    pub key: EncryptionKey,
    /// K in paper
    pub ciphertext: Ciphertext,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    /// k in paper, plaintext of K
    pub plaintext: BigNumber,
    /// rho in paper, nonce of encryption k -> K
    pub nonce: Nonce,
}

// As described in cggmp21 at page 33
/// Prover's first message, obtained by `commit`
#[derive(Debug, Clone)]
pub struct Commitment {
    pub s: BigNumber,
    pub a: BigNumber,
    pub c: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub mu: BigNumber,
    pub r: BigNumber,
    pub gamma: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// `challenge`
pub type Challenge = BigNumber;

// As described in cggmp21 at page 33
/// The ZK proof. Computed by `prove`
#[derive(Debug, Clone)]
pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub z3: BigNumber,
}

pub use crate::common::Aux;

pub mod interactive {
    use crate::unknown_order::BigNumber;
    use rand_core::RngCore;

    use crate::common::{combine, gen_inversible, InvalidProof};

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
    ) -> (Commitment, PrivateCommitment) {
        let two_to_l = BigNumber::from(1) << security.l;
        let two_to_l_plus_e = BigNumber::from(1) << (security.l + security.epsilon);
        let alpha = BigNumber::from_rng(&two_to_l_plus_e, &mut rng);
        let mu = BigNumber::from_rng(&(two_to_l * &aux.rsa_modulo), &mut rng);
        let r = gen_inversible(data.key.n(), &mut rng);
        let gamma = BigNumber::from_rng(&(two_to_l_plus_e * &aux.rsa_modulo), &mut rng);

        let s = combine(&aux.s, &pdata.plaintext, &aux.t, &mu, &aux.rsa_modulo);
        let a = combine(&(data.key.n() + 1), &alpha, &r, data.key.n(), data.key.nn());
        let c = combine(&aux.s, &alpha, &aux.t, &gamma, &aux.rsa_modulo);
        (
            Commitment { s, a, c },
            PrivateCommitment {
                alpha,
                mu,
                r,
                gamma,
            },
        )
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        data: &Data,
        pdata: &PrivateData,
        private_commitment: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Proof {
        let m = crate::unknown_order::Group {
            modulus: data.key.n().clone(),
        };
        let z2 = &m
            * (
                &private_commitment.r,
                &pdata.nonce.modpow(challenge, data.key.n()),
            );
        let z1 = &private_commitment.alpha + (challenge * &pdata.plaintext);
        let z3 = &private_commitment.gamma + (challenge * &private_commitment.mu);
        Proof { z1, z2, z3 }
    }

    /// Verify the proof
    pub fn verify(
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        security: &SecurityParams,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        // check 1
        let pt = &proof.z1 % data.key.n();
        match data.key.encrypt(pt.to_bytes(), Some(proof.z2.clone())) {
            Some((cipher, _nonce)) => {
                if cipher
                    != commitment.a.modmul(
                        &data.ciphertext.modpow(challenge, data.key.nn()),
                        data.key.nn(),
                    )
                {
                    return Err(InvalidProof::EqualityCheckFailed(1));
                }
            }
            None => return Err(InvalidProof::EncryptionFailed),
        }

        let check2 = combine(&aux.s, &proof.z1, &aux.t, &proof.z3, &aux.rsa_modulo)
            == combine(
                &commitment.c,
                &1.into(),
                &commitment.s,
                challenge,
                &aux.rsa_modulo,
            );
        if !check2 {
            return Err(InvalidProof::EqualityCheckFailed(2));
        }

        if proof.z1 > (BigNumber::one() << (security.l + security.epsilon)) {
            return Err(InvalidProof::RangeCheckFailed(3));
        }

        Ok(())
    }

    pub fn challenge<R: RngCore>(security: &SecurityParams, rng: &mut R) -> Challenge {
        BigNumber::from_rng(&security.q, rng)
    }
}

pub mod non_interactive {
    use crate::unknown_order::BigNumber;
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use crate::common::InvalidProof;

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<D, R: RngCore>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        pdata: &PrivateData,
        security: &SecurityParams,
        rng: R,
    ) -> (Commitment, Proof)
    where
        D: Digest<OutputSize = U32>,
    {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng);
        let challenge = challenge(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge);
        (comm, proof)
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<D>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        security: &SecurityParams,
    ) -> Challenge
    where
        D: Digest<OutputSize = U32>,
    {
        use rand_core::SeedableRng;
        let seed = shared_state
            .chain_update(&aux.s.to_bytes())
            .chain_update(&aux.t.to_bytes())
            .chain_update(&aux.rsa_modulo.to_bytes())
            .chain_update(&data.key.to_bytes())
            .chain_update(&data.ciphertext.to_bytes())
            .chain_update(&commitment.s.to_bytes())
            .chain_update(&commitment.a.to_bytes())
            .chain_update(&commitment.c.to_bytes())
            .finalize();
        let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed.into());
        BigNumber::from_rng(&security.q, &mut rng)
    }

    pub fn verify<D>(
        shared_state: D,
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof>
    where
        D: Digest<OutputSize = U32>,
    {
        let challenge = challenge(shared_state, aux, data, commitment, security);
        super::interactive::verify(aux, data, commitment, security, &challenge, proof)
    }
}

#[cfg(test)]
mod test {
    use crate::common::InvalidProof;
    use crate::unknown_order::BigNumber;

    #[test]
    fn passing() {
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
            q: BigNumber::prime(256),
        };
        let private_key = libpaillier::DecryptionKey::random().unwrap();
        let key = libpaillier::EncryptionKey::from(&private_key);
        let plaintext: BigNumber = 228.into();
        let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
        let data = super::Data { key, ciphertext };
        let pdata = super::PrivateData { plaintext, nonce };

        let p = BigNumber::prime(1024);
        let q = BigNumber::prime(1024);
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
            rand_core::OsRng::default(),
        );
        let r = super::non_interactive::verify(
            shared_state,
            &aux,
            &data,
            &commitment,
            &security,
            &proof,
        );
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }
    #[test]
    fn failing() {
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
            q: BigNumber::prime(256),
        };
        let p = BigNumber::prime(1024);
        let q = BigNumber::prime(1024);
        let private_key = libpaillier::DecryptionKey::with_primes(&p, &q).unwrap();
        let key = libpaillier::EncryptionKey::from(&private_key);
        let plaintext: BigNumber = (BigNumber::one() << (security.l + security.epsilon)) + 1;
        let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
        let data = super::Data { key, ciphertext };
        let pdata = super::PrivateData { plaintext, nonce };

        let p = BigNumber::prime(1024);
        let q = BigNumber::prime(1024);
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
            rand_core::OsRng::default(),
        );
        let r = super::non_interactive::verify(
            shared_state,
            &aux,
            &data,
            &commitment,
            &security,
            &proof,
        );
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(InvalidProof::RangeCheckFailed(_)) => (),
            Err(e) => panic!("proof should not fail with {:?}", e),
        }
    }
}
