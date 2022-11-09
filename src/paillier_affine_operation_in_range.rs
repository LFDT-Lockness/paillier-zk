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
//! use paillier_zk::{L, EPSILON};
//! use generic_ec_core::hash_to_curve::Tag;
//! const TAG: Tag = Tag::new_unwrap("application name".as_bytes());
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
//! let ciphertext_mult = g * p::convert_scalar(&plaintext_mult);
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
//! let rng = rand_core::OsRng::default();
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
//! let (commitment, challenge, proof) =
//!     p::compute_proof(TAG, &aux, &data, &pdata, rng).expect("proof failed");
//!
//! // 6. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<C: Curve>(_: &p::Data<C>, _: &p::Commitment<C>, _: &p::Challenge, _: &p::Proof) { todo!() }
//! # fn recv<C: Curve>() -> (p::Data<C>, p::Commitment<C>, p::Challenge, p::Proof) { todo!() }
//! send(&data, &commitment, &challenge, &proof);
//!
//! // 7. Verifier receives the data and the proof and verifies it
//!
//! let (data, commitment, challenge, proof) = recv::<C>();
//! let r = p::verify(&aux, &data, &commitment, &challenge, &proof);
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use crate::unknown_order::BigNumber;
use generic_ec::{Curve, Point, Scalar, hash_to_curve::Tag};
use generic_ec_core::hash_to_curve::HashToCurve;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

use crate::common::{combine, gen_inversible};
pub use crate::common::{InvalidProof, ProtocolError, convert_scalar};
use crate::{EPSILON, L};

/// Public data that both parties know
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
/// Prover's first message, obtained by `commit`
pub struct Commitment<C: Curve> {
    a: BigNumber,
    b_x: Point<C>,
    b_y: BigNumber,
    e: BigNumber,
    s: BigNumber,
    f: BigNumber,
    t: BigNumber,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
pub struct PrivateCommitment {
    alpha: BigNumber,
    beta: BigNumber,
    r: BigNumber,
    r_y: BigNumber,
    gamma: BigNumber,
    m: BigNumber,
    delta: BigNumber,
    mu: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// `challenge`
pub type Challenge = BigNumber;

/// The ZK proof. Computed by `prove`
pub struct Proof {
    z1: BigNumber,
    z2: BigNumber,
    z3: BigNumber,
    z4: BigNumber,
    w: BigNumber,
    w_y: BigNumber,
}

/// Auxiliary data known to both prover and verifier
pub struct Aux {
    /// ring-pedersen parameter
    pub s: BigNumber,
    /// ring-pedersen parameter
    pub t: BigNumber,
    /// N^ in paper
    pub rsa_modulo: BigNumber,
}

/// Create random commitment
pub fn commit<C: Curve, R: RngCore>(
    aux: &Aux,
    data: &Data<C>,
    pdata: &PrivateData,
    mut rng: R,
) -> Result<(Commitment<C>, PrivateCommitment), ProtocolError> {
    let two_to_l = BigNumber::one() << L;
    let two_to_l_e = BigNumber::one() << (L + EPSILON);
    let modulo_l = two_to_l * &aux.rsa_modulo;
    let modulo_l_e = &two_to_l_e * &aux.rsa_modulo;

    let alpha = BigNumber::from_rng(&two_to_l_e, &mut rng);
    let beta = BigNumber::from_rng(&two_to_l_e, &mut rng); // XXX l'
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
        proof.z1 <= &one << (L + EPSILON),
    )?;
    fail_if(
        InvalidProof::RangeCheckFailed(7),
        proof.z2 <= &one << (L + EPSILON), // TODO: L'
    )?;
    Ok(())
}

/// Deterministically compute challenge based on prior known values in protocol
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
            &data.key1.to_bytes(),
            &data.c.to_bytes(),
            &data.d.to_bytes(),
            &data.y.to_bytes(),
            &data.x.to_bytes(true),
            &commitment.a.to_bytes(),
            &commitment.b_x.to_bytes(true),
            &commitment.b_y.to_bytes(),
            &commitment.e.to_bytes(),
            &commitment.s.to_bytes(),
            &commitment.f.to_bytes(),
            &commitment.t.to_bytes(),
        ],
    )
    .map_err(|_| ProtocolError::HashFailed)?;

    Ok(BigNumber::from_slice(scalar.to_be_bytes().as_bytes()))
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
    rng: R,
) -> Result<(Commitment<C>, Challenge, Proof), ProtocolError> {
    let (comm, pcomm) = commit(aux, data, pdata, rng)?;
    let challenge = challenge(tag, aux, data, &comm)?;
    let proof = prove(data, pdata, &pcomm, &challenge);
    Ok((comm, challenge, proof))
}

#[cfg(test)]
mod test {
    use generic_ec::Curve;
    use generic_ec_core::hash_to_curve::HashToCurve;

    use crate::unknown_order::BigNumber;
    use crate::{EPSILON, L};

    fn passing_test<C: Curve + HashToCurve>() {
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);
        let private_key1 = libpaillier::DecryptionKey::random().unwrap();
        let key1 = libpaillier::EncryptionKey::from(&private_key1);
        let g = generic_ec::Point::<C>::generator();
        let plaintext: BigNumber = 228.into();
        let plaintext_orig = BigNumber::from(100);
        let plaintext_mult = BigNumber::from(2);
        let plaintext_add = BigNumber::from(28);
        let (ciphertext, _) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
        let (ciphertext_orig, _) = key0.encrypt(plaintext_orig.to_bytes(), None).unwrap();
        let ciphertext_mult = g * super::convert_scalar(&plaintext_mult);
        let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), None).unwrap();
        let (ciphertext_add_action, nonce) = key0.encrypt(plaintext_add.to_bytes(), None).unwrap();
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

        let p = BigNumber::prime(L + EPSILON + 1);
        let q = BigNumber::prime(L + EPSILON + 1);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) =
            super::compute_proof(tag, &aux, &data, &pdata, rand_core::OsRng::default()).unwrap();
        let r = super::verify(&aux, &data, &commitment, &challenge, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }

    fn failing_test<C: Curve + HashToCurve>() {
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);
        let private_key1 = libpaillier::DecryptionKey::random().unwrap();
        let key1 = libpaillier::EncryptionKey::from(&private_key1);
        let g = generic_ec::Point::<C>::generator();
        let plaintext_orig = BigNumber::from(1337);
        let plaintext_mult = BigNumber::one() << (L + EPSILON) + 1;
        let plaintext_add = BigNumber::one() << (L + EPSILON) + 2;
        let (ciphertext_orig, _) = key0.encrypt(plaintext_orig.to_bytes(), None).unwrap();
        let ciphertext_mult = g * super::convert_scalar(&plaintext_mult);
        let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), None).unwrap();
        let (ciphertext_add_action, nonce) = key0.encrypt(plaintext_add.to_bytes(), None).unwrap();
        // verify that D is obtained from affine transformation of C
        let transformed = key0
            .add(
                &key0.mul(&ciphertext_orig, &plaintext_mult).unwrap(),
                &ciphertext_add_action,
            )
            .unwrap();
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

        let p = BigNumber::prime(L + EPSILON + 1);
        let q = BigNumber::prime(L + EPSILON + 1);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) =
            super::compute_proof(tag, &aux, &data, &pdata, rand_core::OsRng::default()).unwrap();
        let r = super::verify(&aux, &data, &commitment, &challenge, &proof);
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
