#![allow(dead_code)]

use core::ops::Range;
use unknown_order::BigNumber;

trait SigmaProtocol
where
    Self::X: Clone,
    Self::E: Clone,
    Self::Tau: Clone,
{
    type X;
    type Omega;
    type Tau;
    type E;
    type A;
    fn p1(k: Self::X, tau: Self::Tau) -> Self::A;
    fn v1() -> Self::E;
    fn p2(x: Self::X, omega: Self::Omega, tau: Self::Tau, e: Self::E) -> bool;
    fn v2(x: Self::X, a: Self::A, e: Self::E, z: bool) -> bool;
}

fn completeness_prop<P: SigmaProtocol>(x: P::X, omega: P::Omega, tau: P::Tau) -> bool {
    let e = P::v1();
    let a = P::p1(x.clone(), tau.clone());
    let z = P::p2(x.clone(), omega, tau, e.clone());
    P::v2(x, a, e, z) == true
}

fn soundness_prop<P: SigmaProtocol>() {
    // hmmmmm
}

/// Paillier encryption in range
///
/// For public key N, verify that the plaintext of given ciphertext C is in
/// desired range I
fn r_enc(n: BigNumber, i: Range<BigNumber>, c: BigNumber, x: BigNumber, r: BigNumber) -> bool {
    let square_mod = &n * &n;
    let n_plus_1: BigNumber = &n + 1;
    i.contains(&x)
        && c == n_plus_1
            .modpow(&x, &square_mod)
            .modmul(&r.modpow(&n, &square_mod), &square_mod)
}
// One side has 'x' and 'r' that they don't wish to disclose.
// 'C' is the ciphertext of 'x', and is known to all sides.
// This side wants to prove that 'x' is in range 'I'

mod fiat_shamir {
    use rand_core::RngCore;

    pub trait FiatShamir {
        type Data;
        type PrivateData;
        type Commitment;
        type PrivateCommitment;
        type Hash;
        type Proof;
        type Aux;

        fn commit<R: RngCore>(
            aux: &Self::Aux,
            data: &Self::Data,
            private_data: &Self::PrivateData,
            rng: R,
        ) -> (Self::Commitment, Self::PrivateCommitment);
        fn challenge(
            aux: &Self::Aux,
            data: &Self::Data,
            commitment: &Self::Commitment,
        ) -> Self::Hash;
        fn prove(
            aux: &Self::Aux,
            data: &Self::Data,
            private_data: &Self::PrivateData,
            commitment: &Self::Commitment,
            private_commitment: &Self::PrivateCommitment,
            challenge: &Self::Hash,
        ) -> Self::Proof;
        fn verify(
            aux: &Self::Aux,
            data: &Self::Data,
            commitment: &Self::Commitment,
            challenge: &Self::Hash,
            proof: &Self::Proof,
        ) -> Result<(), &'static str>;
    }

    pub fn run_scheme<F: FiatShamir, R: RngCore>(
        aux: &F::Aux,
        data: &F::Data,
        private_data: &F::PrivateData,
        mut rng: R,
    ) -> ((F::Commitment, F::Hash, F::Proof), F::PrivateCommitment) {
        let (commitment, private_commitment) = F::commit(aux, data, private_data, &mut rng);
        let challenge = F::challenge(aux, data, &commitment);
        let proof = F::prove(
            aux,
            data,
            private_data,
            &commitment,
            &private_commitment,
            &challenge,
        );
        ((commitment, challenge, proof), private_commitment)
    }

    pub fn scheme_verify<F>(
        aux: &F::Aux,
        data: &F::Data,
        commitment: &F::Commitment,
        challenge: &F::Hash,
        proof: &F::Proof,
    ) -> Result<(), &'static str>
    where
        F: FiatShamir,
        F::Hash: PartialEq,
    {
        let challenge2 = F::challenge(aux, data, commitment);
        if *challenge != challenge2 {
            Err("challenge doesn't match")
        } else {
            F::verify(aux, data, commitment, challenge, proof)
        }
    }
}

const L: usize = 228;
const EPSILON: usize = 322;

mod paillier_encryption_in_range {
    use crate::{combine, EPSILON, L};
    use rand_core::RngCore;
    use unknown_order::BigNumber;

    pub struct P;

    pub struct Data {
        /// N0 in paper
        key: libpaillier::EncryptionKey,
        /// K in paper
        ciphertext: libpaillier::Ciphertext,
    }
    pub struct PrivateData {
        /// k in paper
        plaintext: BigNumber,
        /// rho in paper
        nonce: libpaillier::Nonce,
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
            let r = crate::gen_inversible(&data.key.n(), &mut rng);
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
}

fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, mut rng: R) -> BigNumber {
    for _ in 0..100 {
        let r = BigNumber::from_rng(modulo, &mut rng);
        if r.gcd(modulo) == 1.into() {
            return r;
        }
    }
    panic!("Attempts exceeded when generating inversible number");
}

mod paillier_affine_operation_in_range {
    use libpaillier::{Ciphertext, EncryptionKey, Nonce};
    use rand_core::RngCore;
    use unknown_order::BigNumber;

    use crate::{combine, gen_inversible, EPSILON, L};

    pub struct P;

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
}

fn combine(
    l: &BigNumber,
    le: &BigNumber,
    r: &BigNumber,
    re: &BigNumber,
    m: &BigNumber,
) -> BigNumber {
    l.modpow(&le, &m).modmul(&r.modpow(&re, &m), &m)
}

fn traverse<'a, I, E>(i: I) -> Result<(), E>
where
    I: IntoIterator<Item = &'a dyn Fn() -> Result<(), E>>,
    E: 'a,
{
    for f in i {
        f()?
    }
    Ok(())
}
