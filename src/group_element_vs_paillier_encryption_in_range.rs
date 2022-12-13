//! ZK-proof, called Пlog* or Rlog* in the CGGMP21 paper.
//!
//! ## Description
//!
//! A party P has a number `X = g ^ x`, with g being a generator of
//! multiplicative group G. P wants to prove to party V that the logarithm of X,
//! i.e. x, is at most L bits.
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
//! use generic_ec_core::hash_to_curve::Tag;
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
//! let power = g * p::convert_scalar(&plaintext);
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 128,
//! };
//!
//! let rng = rand_core::OsRng::default();
//! let data = p::Data { key0, c: ciphertext, x: power, g };
//! let pdata = p::PrivateData { x: plaintext, nonce };
//! let (commitment, challenge, proof) =
//!     p::compute_proof(TAG, &aux, &data, &pdata, &security, rng).expect("proof failed");
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
//! let challenge = p::challenge(TAG, &aux, &data, &commitment).expect("challenge failed");
//! let (data, commitment, proof) = recv::<C>();
//! p::verify(&aux, &data, &commitment, &security, &challenge, &proof);
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use crate::{
    common::{combine, gen_inversible, InvalidProof, ProtocolError},
    unknown_order::BigNumber,
};
use generic_ec::{hash_to_curve::Tag, Curve, Point, Scalar};
use generic_ec_core::hash_to_curve::HashToCurve;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

pub use crate::common::convert_scalar;

pub struct SecurityParams {
    /// l in paper, bit size of plaintext
    pub l: usize,
    /// Epsilon in paper, range extension and security parameter for x
    pub epsilon: usize,
}

pub struct Data<C: Curve> {
    pub key0: EncryptionKey,
    pub c: Ciphertext,
    pub x: Point<C>,
    /// A generator in group
    pub g: Point<C>,
}

pub struct PrivateData {
    pub x: BigNumber,
    pub nonce: Nonce,
}

pub struct Commitment<C: Curve> {
    pub s: BigNumber,
    pub a: Ciphertext,
    pub y: Point<C>,
    pub d: BigNumber,
}

pub struct PrivateCommitment {
    pub alpha: BigNumber,
    pub mu: BigNumber,
    pub r: Nonce,
    pub gamma: BigNumber,
}

pub type Challenge = BigNumber;

pub struct Proof {
    pub z1: BigNumber,
    pub z2: BigNumber,
    pub z3: BigNumber,
}

pub use crate::common::Aux;

/// Create random commitment
pub fn commit<C: Curve, R: RngCore>(
    aux: &Aux,
    data: &Data<C>,
    pdata: &PrivateData,
    security: &SecurityParams,
    mut rng: R,
) -> Result<(Commitment<C>, PrivateCommitment), ProtocolError> {
    let two_to_l = BigNumber::one() << security.l;
    let two_to_l_e = BigNumber::one() << (security.l + security.epsilon);
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

pub fn challenge<C: Curve + HashToCurve>(
    tag: Tag,
    aux: &Aux,
    data: &Data<C>,
    commitment: &Commitment<C>,
) -> Result<Challenge, ProtocolError> {
    use generic_ec::hash_to_curve::FromHash;
    let scalar = Scalar::<C>::hash_concat(
        tag,
        &[
            aux.s.to_bytes().as_ref(), // hint for array to become [&[u8]]
            &aux.t.to_bytes(),
            &aux.rsa_modulo.to_bytes(),
            &data.key0.to_bytes(),
            &data.c.to_bytes(),
            &data.x.to_bytes(true),
            &commitment.s.to_bytes(),
            &commitment.a.to_bytes(),
            &commitment.y.to_bytes(true),
            &commitment.d.to_bytes(),
        ],
    )
    .map_err(|_| ProtocolError::HashFailed)?;

    Ok(BigNumber::from_slice(scalar.to_be_bytes().as_bytes()))
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
        proof.z1 <= one << (security.l + security.epsilon),
        InvalidProof::RangeCheckFailed(4),
    )?;

    Ok(())
}

/// Compute proof for the given data, producing random commitment and
/// deriving determenistic challenge.
///
/// Obtained from the above interactive proof via Fiat-Shamir heuristic.
pub fn compute_proof<C: Curve + HashToCurve, R: RngCore>(
    tag: Tag,
    aux: &Aux,
    data: &Data<C>,
    pdata: &PrivateData,
    security: &SecurityParams,
    rng: R,
) -> Result<(Commitment<C>, Challenge, Proof), ProtocolError> {
    let (comm, pcomm) = commit(aux, data, pdata, security, rng)?;
    let challenge = challenge(tag, aux, data, &comm)?;
    let proof = prove(data, pdata, &pcomm, &challenge);
    Ok((comm, challenge, proof))
}

#[cfg(test)]
mod test {
    use generic_ec::Curve;
    use generic_ec_core::hash_to_curve::HashToCurve;
    use libpaillier::unknown_order::BigNumber;

    use crate::common::convert_scalar;

    fn passing_test<C: Curve + HashToCurve>() {
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
        };
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(228);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
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

        let p = BigNumber::prime(1024);
        let q = BigNumber::prime(1024);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());

        let (commitment, challenge, proof) =
            super::compute_proof(tag, &aux, &data, &pdata, &security, rand_core::OsRng::default()).unwrap();
        let r = super::verify(&aux, &data, &commitment, &security, &challenge, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }

    fn failing_test<C: Curve + HashToCurve>() {
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 128,
        };
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(1) << (security.l + security.epsilon + 1);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
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

        let p = BigNumber::prime(1024);
        let q = BigNumber::prime(1024);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());

        let (commitment, challenge, proof) =
            super::compute_proof(tag, &aux, &data, &pdata, &security, rand_core::OsRng::default()).unwrap();
        let r = super::verify(&aux, &data, &commitment, &security, &challenge, &proof);
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(_) => (),
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
}
