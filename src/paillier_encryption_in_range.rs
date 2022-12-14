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
//! const TAG: Tag = Tag::new_unwrap("application name".as_bytes());
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
//! let (commitment, challenge, proof) =
//!     p::compute_proof(TAG, &aux, &data, &pdata, &security, rng);
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Proof) { todo!() }
//! # fn recv() -> (p::Data, p::Commitment, p::Proof) { todo!() }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! let challenge = p::challenge(TAG, &aux, &data, &commitment, &security);
//! let (data, commitment, proof) = recv();
//! p::verify(&aux, &data, &commitment, &security, &challenge, &proof);
//!
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use crate::unknown_order::BigNumber;
use generic_ec::hash_to_curve::Tag;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

pub use crate::common::InvalidProof;
use crate::common::{combine, gen_inversible};

pub struct SecurityParams {
    /// l in paper, bit size of plaintext
    pub l: usize,
    /// Epsilon in paper, range extension and security parameter for x
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: BigNumber,
}

/// Public data that both parties know
pub struct Data {
    /// N0 in paper, public key that k -> K was encrypted on
    pub key: EncryptionKey,
    /// K in paper
    pub ciphertext: Ciphertext,
}

/// Private data of prover
pub struct PrivateData {
    /// k in paper, plaintext of K
    pub plaintext: BigNumber,
    /// rho in paper, nonce of encryption k -> K
    pub nonce: Nonce,
}

// As described in cggmp21 at page 33
/// Prover's first message, obtained by `commit`
pub struct Commitment {
    pub s: BigNumber,
    pub a: BigNumber,
    pub c: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
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
pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub z3: BigNumber,
}

pub use crate::common::Aux;

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
#[allow(clippy::just_underscores_and_digits)]
fn prove(
    data: &Data,
    pdata: &PrivateData,
    private_commitment: &PrivateCommitment,
    challenge: &Challenge,
) -> Proof {
    let m = crate::unknown_order::Group {
        modulus: data.key.n().clone(),
    };
    let _2 = &m
        * (
            &private_commitment.r,
            &pdata.nonce.modpow(challenge, data.key.n()),
        );
    let _1 = &private_commitment.alpha + (challenge * &pdata.plaintext);
    let _3 = &private_commitment.gamma + (challenge * &private_commitment.mu);
    Proof {
        z1: _1,
        z2: _2,
        z3: _3,
    }
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

/// Deterministically compute challenge based on prior known values in protocol
pub fn challenge(
    tag: Tag,
    aux: &Aux,
    data: &Data,
    commitment: &Commitment,
    security: &SecurityParams,
) -> Challenge {
    crate::common::hash2field::hash_to_field(
        tag,
        &security.q,
        &[
            &aux.s.to_bytes(),
            &aux.t.to_bytes(),
            &aux.rsa_modulo.to_bytes(),
            &data.key.to_bytes(),
            &data.ciphertext.to_bytes(),
            &commitment.s.to_bytes(),
            &commitment.a.to_bytes(),
            &commitment.c.to_bytes(),
        ],
    )
}

/// Compute proof for the given data, producing random commitment and
/// deriving determenistic challenge.
///
/// Obtained from the above interactive proof via Fiat-Shamir heuristic.
pub fn compute_proof<R: RngCore>(
    tag: Tag,
    aux: &Aux,
    data: &Data,
    pdata: &PrivateData,
    security: &SecurityParams,
    rng: R,
) -> (Commitment, Challenge, Proof) {
    let (comm, pcomm) = commit(aux, data, pdata, security, rng);
    let challenge = challenge(tag, aux, data, &comm, security);
    let proof = prove(data, pdata, &pcomm, &challenge);
    (comm, challenge, proof)
}

#[cfg(test)]
mod test {
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

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) = super::compute_proof(
            tag,
            &aux,
            &data,
            &pdata,
            &security,
            rand_core::OsRng::default(),
        );
        let r = super::verify(&aux, &data, &commitment, &security, &challenge, &proof);
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

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) = super::compute_proof(
            tag,
            &aux,
            &data,
            &pdata,
            &security,
            rand_core::OsRng::default(),
        );
        let r = super::verify(&aux, &data, &commitment, &security, &challenge, &proof);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(_) => (),
        }
    }
}
