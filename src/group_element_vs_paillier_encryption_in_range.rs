use crate::{unknown_order::BigNumber, common::{combine, convert_scalar, gen_inversible, ProtocolError, InvalidProof}, EPSILON, L};
use generic_ec::{Curve, Point, hash_to_curve::Tag, Scalar};
use generic_ec_core::hash_to_curve::HashToCurve;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

pub struct Data<C: Curve> {
    pub key0: EncryptionKey,
    pub c: Ciphertext,
    pub x: Point<C>,
}

pub struct PrivateData {
    pub x: BigNumber,
    pub nonce: Nonce,
}

pub struct Commitment<C: Curve> {
    s: BigNumber,
    a: Ciphertext,
    y: Point<C>,
    d: BigNumber,
}

pub struct PrivateCommitment {
    alpha: BigNumber,
    mu: BigNumber,
    r: Nonce,
    gamma: BigNumber,
}

pub type Challenge = BigNumber;

pub struct Proof {
    z1: BigNumber,
    z2: BigNumber,
    z3: BigNumber,
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
    let mu = BigNumber::from_rng(&modulo_l, &mut rng);
    let r = gen_inversible(data.key0.n(), &mut rng);
    let gamma = BigNumber::from_rng(&modulo_l_e, &mut rng);

    let (a, _) = data.key0.encrypt(alpha.to_bytes(), Some(r.clone())).ok_or(ProtocolError::EncryptionFailed)?;

    let commitment = Commitment {
        s: combine(&aux.s, &pdata.x, &aux.t, &mu, &aux.rsa_modulo),
        a,
        y: Point::<C>::generator() * convert_scalar(&alpha),
        d: combine(&aux.s, &alpha, &aux.t, &gamma, &aux.rsa_modulo),
    };
    let private_commitment = PrivateCommitment {
        alpha, mu, r, gamma
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
        z2: combine(&pcomm.r, &BigNumber::one(), &pdata.nonce, &challenge, data.key0.n()),
        z3: &pcomm.gamma + challenge * &pcomm.mu,
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
    fn fail_if(b: bool, msg: InvalidProof) -> Result<(), InvalidProof> {
        if b {
            Ok(())
        } else {
            Err(msg)
        }
    }
    // Three equality checks and one range check
    {
        let (lhs, _) = data.key0.encrypt(proof.z1.to_bytes(), Some(proof.z2.clone())).ok_or(InvalidProof::EncryptionFailed)?;
        let rhs = combine(&commitment.a, &one, &data.c, &challenge, &data.key0.nn());
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(1))?;
    }
    {
        let lhs = Point::<C>::generator() * convert_scalar(&proof.z1);
        let rhs = commitment.y + data.x * convert_scalar(challenge);
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(2))?;
    }
    {
        let lhs = combine(&aux.s, &proof.z1, &aux.t, &proof.z3, &aux.rsa_modulo);
        let rhs = combine(&commitment.d, &one, &commitment.s, &challenge, &aux.rsa_modulo);
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(3))?;
    }
    fail_if( proof.z1 <= one << (L + EPSILON), InvalidProof::RangeCheckFailed(4) )?;

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
    use libpaillier::unknown_order::BigNumber;

    use crate::{common::convert_scalar, L, EPSILON};

    fn passing_test<C: Curve + HashToCurve>() {
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(228);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
        let x = generic_ec::Point::<C>::generator() * convert_scalar(&plaintext);

        let data = super::Data {
            key0,
            c: ciphertext,
            x,
        };
        let pdata = super::PrivateData {
            x: plaintext,
            nonce,
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

        let plaintext = BigNumber::from(1) << ( L + EPSILON ) + 1;
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
        let x = generic_ec::Point::<C>::generator() * convert_scalar(&plaintext);

        let data = super::Data {
            key0,
            c: ciphertext,
            x,
        };
        let pdata = super::PrivateData {
            x: plaintext,
            nonce,
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
