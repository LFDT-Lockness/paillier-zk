//! Paillier encryption in range, Пenc
//!
//! For public key N, verify that the plaintext of given ciphertext C is in
//! desired range I.  
//! Range is given as bitsize by parameters epsilon and l

use libpaillier::{EncryptionKey, Ciphertext, Nonce};
use rand_core::RngCore;
use unknown_order::BigNumber;

use crate::{combine, EPSILON, L, gen_inversible};

pub struct P;

pub struct Data {
    /// N0 in paper, public key that k -> K was encrypted on
    key: EncryptionKey,
    /// K in paper
    ciphertext: Ciphertext,
}
pub struct PrivateData {
    /// k in paper, plaintext of K
    plaintext: BigNumber,
    /// rho in paper, nonce of encryption k -> K
    nonce: Nonce,
}
/// As described in cggmp21 at page 33
pub struct Commitment {
    s: BigNumber,
    a: BigNumber,
    c: BigNumber,
}
pub struct PrivateCommitment {
    alpha: BigNumber,
    mu: BigNumber,
    r: BigNumber,
    gamma: BigNumber,
}
pub type Hash = BigNumber;
/// As described in cggmp21 at page 33
pub struct Proof {
    _1: BigNumber,
    _2: BigNumber,
    _3: BigNumber,
}
pub struct Aux {
    /// ring-pedersen parameter
    s: BigNumber,
    /// ring-pedersen parameter
    t: BigNumber,
    /// N^ in paper
    rsa_modulo: BigNumber,
}

impl crate::fiat_shamir::FiatShamir for P {
    type Data = Data;
    type PrivateData = PrivateData;
    type Commitment = Commitment;
    type PrivateCommitment = PrivateCommitment;
    type Hash = Hash;
    type Proof = Proof;
    type Aux = Aux;

    fn commit<R: RngCore>(
        aux: &Aux,
        data: &Data,
        pdata: &PrivateData,
        mut rng: R,
    ) -> (Commitment, PrivateCommitment) {
        let two_to_l = BigNumber::from(1) << L;
        let two_to_l_plus_e = BigNumber::from(1) << (L + EPSILON);
        let alpha = BigNumber::from_rng(&two_to_l_plus_e, &mut rng);
        let mu = BigNumber::from_rng(&(two_to_l * &aux.rsa_modulo), &mut rng);
        let r = gen_inversible(&data.key.n(), &mut rng);
        let gamma = BigNumber::from_rng(&(two_to_l_plus_e * &aux.rsa_modulo), &mut rng);

        let s = combine(&aux.s, &pdata.plaintext, &aux.t, &mu, &aux.rsa_modulo);
        let a = combine(
            &(data.key.n() + 1),
            &alpha,
            &r,
            &data.key.n(),
            &data.key.nn(),
        );
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

    fn prove(
        _: &Aux,
        data: &Data,
        pdata: &PrivateData,
        _: &Commitment,
        private_commitment: &PrivateCommitment,
        challenge: &Hash,
    ) -> Self::Proof {
        let m = unknown_order::Group {
            modulus: data.key.n().clone(),
        };
        let _2 = &m
            * (
                &private_commitment.r,
                &pdata.nonce.modpow(challenge, &data.key.n()),
            );
        let _1 = &private_commitment.alpha + (challenge * &pdata.plaintext);
        let _3 = &private_commitment.gamma + (challenge * &private_commitment.mu);
        Proof { _1, _2, _3 }
    }

    fn verify(
        aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        challenge: &Hash,
        proof: &Proof,
    ) -> Result<(), &'static str> {
        let check1 = || {
            let pt = &proof._1 % data.key.n();
            match data.key.encrypt(pt.to_bytes(), Some(proof._2.clone())) {
                Some((cipher, _nonce)) => {
                    if cipher
                        != commitment.a.modmul(
                            &data.ciphertext.modpow(challenge, data.key.nn()),
                            data.key.nn(),
                        )
                    {
                        Err("check1 failed")
                    } else {
                        Ok(())
                    }
                }
                None => Err("encrypt failed"),
            }
        };
        fn fail_if(b: bool, msg: &'static str) -> Result<(), &'static str> {
            if b {
                Ok(())
            } else {
                Err(msg)
            }
        }
        let check2 = || {
            fail_if(
                combine(&aux.s, &proof._1, &aux.t, &proof._3, &aux.rsa_modulo)
                    == combine(
                        &commitment.c,
                        &1.into(),
                        &commitment.s,
                        challenge,
                        &aux.rsa_modulo,
                    ),
                "check2",
            )
        };
        let check3 = || fail_if(proof._1 <= (BigNumber::one() << (L + EPSILON)), "check3");

        // need to explicitly erase type
        let c1: &dyn Fn() -> Result<(), &'static str> = &check1;
        crate::traverse([c1, &check2, &check3])
    }

    fn challenge(aux: &Aux, data: &Data, commitment: &Commitment) -> Hash {
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

        BigNumber::from_slice(digest.finalize())
    }
}

#[cfg(test)]
mod test {
    use super::{EPSILON, L};
    use unknown_order::BigNumber;

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

        let ((commitment, challenge, proof), _pcomm) =
            crate::fiat_shamir::run_scheme::<super::P, _>(
                &aux,
                &data,
                &pdata,
                rand_core::OsRng::default(),
            );
        let r = crate::fiat_shamir::scheme_verify::<super::P>(
            &aux,
            &data,
            &commitment,
            &challenge,
            &proof,
        );
        match r {
            Ok(()) => (),
            Err(e) => panic!("{}", e),
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

        let ((commitment, challenge, proof), _pcomm) =
            crate::fiat_shamir::run_scheme::<super::P, _>(
                &aux,
                &data,
                &pdata,
                rand_core::OsRng::default(),
            );
        let r = crate::fiat_shamir::scheme_verify::<super::P>(
            &aux,
            &data,
            &commitment,
            &challenge,
            &proof,
        );
        match r {
            Ok(()) => panic!("proof should not pass"),
            Err(_) => (),
        }
    }
}

