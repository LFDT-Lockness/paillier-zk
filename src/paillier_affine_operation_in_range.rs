//! Paillier operation with group commitment in range, Пaff-g
//!
//! Given a lot of ciphertexts and keys, having C = kx + c, verify that k and c
//! are in desired ranges I and J.  
//! Ranges are given as bitsize by parameters epsilon, l, and l'

use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;
use unknown_order::BigNumber;

use crate::{combine, gen_inversible, EPSILON, L};

pub struct P;

/// Parameterized by G, which is the group that hides the x parameter. I think
/// this can only work with a group being a group of congruence classes in Z,
/// reinterpreted to Z when needed. Certainly this is the only way we will use
/// it, so I will remove this parameter later.
pub struct Data<G> {
    /// Group generator
    g: G,
    /// Group rank
    q: BigNumber,
    /// N0 in paper, public key that C was encrypted on
    key0: EncryptionKey,
    /// N1 in paper, public key that y -> Y was encrypted on
    key1: EncryptionKey,
    /// C or C0 in paper, some data encrypted on N0
    c: Ciphertext,
    /// D or C in paper, result of affine transformation of C0 with x and y
    d: BigNumber,
    /// Y in paper, y encrypted on N1
    y: Ciphertext,
    /// X in paper, obtained as g^x
    x: BigNumber,
}
pub struct PrivateData {
    /// x or epsilon in paper, preimage of X
    x: BigNumber,
    /// y or delta in paper, preimage of Y
    y: BigNumber,
    /// rho in paper, nonce in encryption of y for additive action
    nonce: Nonce,
    /// rho_y in paper, nonce in encryption of y to obtain Y
    nonce_y: Nonce,
}
/// As described in cggmp21 at page 35
pub struct Commitment<G> {
    a: BigNumber,
    b_x: G,
    b_y: BigNumber,
    e: BigNumber,
    s: BigNumber,
    f: BigNumber,
    t: BigNumber,
}
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
pub type Hash = BigNumber;
pub struct Proof {
    z1: BigNumber,
    z2: BigNumber,
    z3: BigNumber,
    z4: BigNumber,
    w: BigNumber,
    w_y: BigNumber,
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
    // It seems this is impossible to write with anything other than
    // BigNumber as generic parameter
    type Data = Data<BigNumber>;
    type PrivateData = PrivateData;
    type Commitment = Commitment<BigNumber>;
    type PrivateCommitment = PrivateCommitment;
    type Hash = Hash;
    type Proof = Proof;
    type Aux = Aux;

    fn commit<R: RngCore>(
        aux: &Aux,
        data: &Self::Data,
        pdata: &PrivateData,
        mut rng: R,
    ) -> (Self::Commitment, PrivateCommitment) {
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

        let a_add = data.key0.encrypt(beta.to_bytes(), Some(r.clone())).unwrap().0;
        let c_to_alpha = data.key0.mul(&data.c, &alpha).unwrap();
        let a = data.key0.add(&c_to_alpha, &a_add).unwrap();
        let commitment = Commitment {
            a,
            b_x: data.g.modpow(&alpha, &data.q),
            b_y: data.key1.encrypt(beta.to_bytes(), Some(r_y.clone())).unwrap().0,
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
        (commitment, private_commitment)
    }

    fn prove(
        _: &Aux,
        data: &Self::Data,
        pdata: &PrivateData,
        _: &Self::Commitment,
        pcomm: &PrivateCommitment,
        challenge: &Hash,
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
                &challenge,
                &data.key0.n(),
            ),
            w_y: combine(
                &pcomm.r_y,
                &BigNumber::one(),
                &pdata.nonce_y,
                &challenge,
                &data.key1.n(),
            ),
        }
    }

    fn verify(
        aux: &Aux,
        data: &Self::Data,
        commitment: &Self::Commitment,
        challenge: &Hash,
        proof: &Proof,
    ) -> Result<(), &'static str> {
        let one = BigNumber::one();
        fn fail_if(msg: &'static str, b: bool) -> Result<(), &'static str> {
            if b {
                Ok(())
            } else {
                Err(msg)
            }
        }
        let check1 = || {
            let enc = data
                .key0
                .encrypt(proof.z2.to_bytes(), Some(proof.w.clone()))
                .unwrap()
                .0;
            let lhs = data
                .key0
                .add(&data.key0.mul(&data.c, &proof.z1).unwrap(), &enc)
                .unwrap();
            let rhs = combine(&commitment.a, &one, &data.d, &challenge, &data.key0.nn());
            fail_if("check1", lhs == rhs)
        };
        let check2 = || {
            let lhs = data.g.modpow(&proof.z1, &data.q);
            let rhs = combine(&commitment.b_x, &one, &data.x, &challenge, &data.q);
            fail_if("check2", lhs == rhs)
        };
        let check3 = || {
            let lhs = data
                .key1
                .encrypt(proof.z2.to_bytes(), Some(proof.w_y.clone()))
                .unwrap()
                .0;
            let rhs = combine(&commitment.b_y, &one, &data.y, &challenge, &data.key1.nn());
            fail_if("check3", lhs == rhs)
        };
        let check4 = || {
            fail_if(
                "check4",
                combine(&aux.s, &proof.z1, &aux.t, &proof.z3, &aux.rsa_modulo)
                    == combine(
                        &commitment.e,
                        &one,
                        &commitment.s,
                        &challenge,
                        &aux.rsa_modulo,
                    ),
            )
        };
        let check5 = || {
            fail_if(
                "check5",
                combine(&aux.s, &proof.z2, &aux.t, &proof.z4, &aux.rsa_modulo)
                    == combine(
                        &commitment.f,
                        &one,
                        &commitment.t,
                        &challenge,
                        &aux.rsa_modulo,
                    ),
            )
        };
        let check6 = || {
            fail_if(
                "range check6",
                proof.z1 <= &one << (L + EPSILON),
            )
        };
        let check7 = || {
            fail_if(
                "range check7",
                proof.z2 <= &one << (L + EPSILON), // TODO: L'
            )
        };
        let c1: &dyn Fn() -> Result<(), &'static str> = &check1;
        crate::traverse([c1, &check2, &check3, &check4, &check5, &check6, &check7])
    }

    fn challenge(
        _aux: &Aux,
        _data: &Data<BigNumber>,
        _commitment: &Commitment<BigNumber>,
    ) -> Hash {
        1.into()
    }
}

#[cfg(test)]
mod test {
    use unknown_order::BigNumber;

    use crate::{EPSILON, L};

    #[test]
    fn passing() {
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);
        let private_key1 = libpaillier::DecryptionKey::random().unwrap();
        let key1 = libpaillier::EncryptionKey::from(&private_key1);
        let plaintext: BigNumber = 228.into();
        let plaintext_orig = BigNumber::from(100);
        let plaintext_mult = BigNumber::from(2);
        let plaintext_add = BigNumber::from(28);
        let q = BigNumber::from(1_000_000_007);
        let g = BigNumber::from(2);
        // verify that g is generator in Z/q
        assert_eq!(g.gcd(&q), 1.into());
        let (ciphertext, _) = key0.encrypt(plaintext.to_bytes(), None).unwrap();
        let (ciphertext_orig, _) =
            key0.encrypt(plaintext_orig.to_bytes(), None).unwrap();
        let ciphertext_mult = g.modpow(&plaintext_mult, &q);
        let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), None).unwrap();
        let (ciphertext_add_action, nonce) =
            key0.encrypt(plaintext_add.to_bytes(), None).unwrap();
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
            g,
            q,
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
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);
        let private_key1 = libpaillier::DecryptionKey::random().unwrap();
        let key1 = libpaillier::EncryptionKey::from(&private_key1);
        let plaintext_orig = BigNumber::from(1337);
        let plaintext_mult = BigNumber::one() << (L + EPSILON) + 1;
        let plaintext_add = BigNumber::one() << (L + EPSILON) + 2;
        let q = BigNumber::from(1_000_000_007);
        let g = BigNumber::from(2);
        // verify that g is generator in Z/q
        assert_eq!(g.gcd(&q), 1.into());
        let (ciphertext_orig, _) =
            key0.encrypt(plaintext_orig.to_bytes(), None).unwrap();
        let ciphertext_mult = g.modpow(&plaintext_mult, &q);
        let (ciphertext_add, nonce_y) = key1.encrypt(plaintext_add.to_bytes(), None).unwrap();
        let (ciphertext_add_action, nonce) =
            key0.encrypt(plaintext_add.to_bytes(), None).unwrap();
        // verify that D is obtained from affine transformation of C
        let transformed = key0
            .add(
                &key0.mul(&ciphertext_orig, &plaintext_mult).unwrap(),
                &ciphertext_add_action,
            )
            .unwrap();
        let data = super::Data {
            g,
            q,
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

