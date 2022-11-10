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
//! use paillier_zk::{L, EPSILON};
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
//!
//! let p = BigNumber::prime(L + EPSILON + 1);
//! let q = BigNumber::prime(L + EPSILON + 1);
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
//! let key = libpaillier::EncryptionKey::from(&private_key);
//!
//! // 2. Setup: prover has some plaintext and encrypts it
//!
//! let plaintext: BigNumber = 228.into();
//! let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most `L` bits:
//!
//! let rng = rand_core::OsRng::default();
//! let data = p::Data { key, ciphertext };
//! let pdata = p::PrivateData { plaintext, nonce };
//! let (commitment, challenge, proof) =
//!     p::compute_proof(&aux, &data, &pdata, rng);
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Challenge, _: &p::Proof) { todo!() }
//! # fn recv() -> (p::Data, p::Commitment, p::Challenge, p::Proof) { todo!() }
//! send(&data, &commitment, &challenge, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! let (data, commitment, challenge, proof) = recv();
//! p::verify(&aux, &data, &commitment, &challenge, &proof);
//!
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use crate::unknown_order::BigNumber;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

use crate::common::{combine, gen_inversible};
use crate::{EPSILON, L};
pub use crate::common::InvalidProof;

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
    s: BigNumber,
    a: BigNumber,
    c: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
pub struct PrivateCommitment {
    alpha: BigNumber,
    mu: BigNumber,
    r: BigNumber,
    gamma: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// `challenge`
pub type Challenge = BigNumber;

// As described in cggmp21 at page 33
/// The ZK proof. Computed by `prove`
pub struct Proof {
    _1: BigNumber,
    _2: BigNumber,
    _3: BigNumber,
}

pub use crate::common::Aux;

/// Create random commitment
pub fn commit<R: RngCore>(
    aux: &Aux,
    data: &Data,
    pdata: &PrivateData,
    mut rng: R,
) -> (Commitment, PrivateCommitment) {
    let two_to_l = BigNumber::from(1) << L;
    let two_to_l_plus_e = BigNumber::from(1) << (L + EPSILON);
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
    Proof { _1, _2, _3 }
}

/// Verify the proof
pub fn verify(
    aux: &Aux,
    data: &Data,
    commitment: &Commitment,
    challenge: &Challenge,
    proof: &Proof,
) -> Result<(), InvalidProof> {
    // check 1
    let pt = &proof._1 % data.key.n();
    match data.key.encrypt(pt.to_bytes(), Some(proof._2.clone())) {
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

    let check2 = combine(&aux.s, &proof._1, &aux.t, &proof._3, &aux.rsa_modulo)
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

    if proof._1 > (BigNumber::one() << (L + EPSILON)) {
        return Err(InvalidProof::RangeCheckFailed(3));
    }

    Ok(())
}

/// Deterministically compute challenge based on prior known values in protocol
pub fn challenge(aux: &Aux, data: &Data, commitment: &Commitment) -> Challenge {
    use sha2::Digest;
    let mut digest = sha2::Sha512::new();

    digest.update(aux.s.to_bytes());
    digest.update(aux.t.to_bytes());
    digest.update(aux.rsa_modulo.to_bytes());

    digest.update(data.key.to_bytes());
    digest.update(data.ciphertext.to_bytes());

    digest.update(commitment.s.to_bytes());
    digest.update(commitment.a.to_bytes());
    digest.update(commitment.c.to_bytes());

    // FIXME: hash to bignumber
    BigNumber::from_slice(digest.finalize())
}

/// Compute proof for the given data, producing random commitment and
/// deriving determenistic challenge.
///
/// Obtained from the above interactive proof via Fiat-Shamir heuristic.
pub fn compute_proof<R: RngCore>(
    aux: &Aux,
    data: &Data,
    pdata: &PrivateData,
    rng: R,
) -> (Commitment, Challenge, Proof) {
    let (comm, pcomm) = commit(aux, data, pdata, rng);
    let challenge = challenge(aux, data, &comm);
    let proof = prove(data, pdata, &pcomm, &challenge);
    (comm, challenge, proof)
}

#[cfg(test)]
mod test {
    use super::{EPSILON, L};
    use crate::unknown_order::BigNumber;

    #[test]
    fn passing() {
        let private_key = libpaillier::DecryptionKey::random().unwrap();
        let key = libpaillier::EncryptionKey::from(&private_key);
        let plaintext: BigNumber = 228.into();
        let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
        let data = super::Data { key, ciphertext };
        let pdata = super::PrivateData { plaintext, nonce };

        let p = BigNumber::prime(L + EPSILON + 1);
        let q = BigNumber::prime(L + EPSILON + 1);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let (commitment, challenge, proof) =
            super::compute_proof(&aux, &data, &pdata, rand_core::OsRng::default());
        let r = super::verify(&aux, &data, &commitment, &challenge, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }
    #[test]
    fn failing() {
        let p = BigNumber::prime(L + 1);
        let q = BigNumber::prime(EPSILON + 1);
        let private_key = libpaillier::DecryptionKey::with_primes(&p, &q).unwrap();
        let key = libpaillier::EncryptionKey::from(&private_key);
        let plaintext: BigNumber = (BigNumber::one() << (L + EPSILON)) + 1;
        let (ciphertext, nonce) = key.encrypt(plaintext.to_bytes(), None).unwrap();
        let data = super::Data { key, ciphertext };
        let pdata = super::PrivateData { plaintext, nonce };

        let p = BigNumber::prime(L + EPSILON + 1);
        let q = BigNumber::prime(L + EPSILON + 1);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let (commitment, challenge, proof) =
            super::compute_proof(&aux, &data, &pdata, rand_core::OsRng::default());
        let r = super::verify(&aux, &data, &commitment, &challenge, &proof);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(_) => (),
        }
    }
}
