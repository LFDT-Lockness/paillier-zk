#![allow(dead_code)]

use unknown_order::BigNumber;
use core::ops::Range;

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
    i.contains(&x) &&
        c == n_plus_1.modpow(&x, &square_mod).modmul(&r.modpow(&n, &square_mod), &square_mod)
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

        fn commit<R: RngCore>(aux: &Self::Aux, data: &Self::Data, private_data: &Self::PrivateData, rng: R) -> (Self::Commitment, Self::PrivateCommitment);
        fn challenge(aux: &Self::Aux, data: &Self::Data, commitment: &Self::Commitment) -> Self::Hash;
        fn prove(aux: &Self::Aux, data: &Self::Data, private_data: &Self::PrivateData, commitment: &Self::Commitment, private_commitment: &Self::PrivateCommitment, challenge: &Self::Hash) -> Self::Proof;
        fn verify(aux: &Self::Aux, data: &Self::Data, commitment: &Self::Commitment, challenge: &Self::Hash, proof: &Self::Proof) -> Result<(), &'static str>;
    }

    pub fn run_scheme<F: FiatShamir, R: RngCore>(aux: &F::Aux, data: &F::Data, private_data: &F::PrivateData, mut rng: R) -> ((F::Commitment, F::Hash, F::Proof), F::PrivateCommitment) {
        let (commitment, private_commitment) = F::commit(aux, data, private_data, &mut rng);
        let challenge = F::challenge(aux, data, &commitment);
        let proof = F::prove(aux, data, private_data, &commitment, &private_commitment, &challenge);
        ((commitment, challenge, proof), private_commitment)
    }

    pub fn scheme_verify<F>(aux: &F::Aux, data: &F::Data, commitment: &F::Commitment, challenge: &F::Hash, proof: &F::Proof) -> Result<(), &'static str>
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

mod paillier_encryption_in_range {
    use rand_core::RngCore;
    use unknown_order::BigNumber;

    struct P;

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

    const L: usize = 228;
    const EPSILON: usize = 322;

    impl crate::fiat_shamir::FiatShamir for P {
        type Data = Data;
        type PrivateData = PrivateData;
        type Commitment = Commitment;
        type PrivateCommitment = PrivateCommitment;
        type Hash = Hash;
        type Proof = Proof;
        type Aux = Aux;

        fn commit<R: RngCore>(aux: &Aux, data: &Data, pdata: &PrivateData, mut rng: R) -> (Self::Commitment, Self::PrivateCommitment) {
            let two_to_l = BigNumber::from(1) << L;
            let two_to_l_plus_e = BigNumber::from(1) << (L + EPSILON);
            let alpha = BigNumber::from_rng(&two_to_l_plus_e, &mut rng);
            let mu = BigNumber::from_rng(&(two_to_l * &aux.rsa_modulo), &mut rng);
            let r = crate::gen_inversible(&data.key.n(), &mut rng);
            let gamma = BigNumber::from_rng(&(two_to_l_plus_e * &aux.rsa_modulo), &mut rng);

            let s = combine(&aux.s, &pdata.plaintext, &aux.t, &mu, &aux.rsa_modulo);
            let a = combine(&(data.key.n() + 1), &alpha, &r, &data.key.n(), &data.key.nn());
            let c = combine(&aux.s, &alpha, &aux.t, &gamma, &aux.rsa_modulo);
            (Commitment {s, a, c}, PrivateCommitment {alpha, mu, r, gamma})
        }

        fn prove(_: &Aux, data: &Data, pdata: &PrivateData, _: &Commitment, private_commitment: &PrivateCommitment, challenge: &Hash) -> Self::Proof {
            let m = unknown_order::Group { modulus: data.key.n().clone() };
            let _2 = &m * (&private_commitment.r, &pdata.nonce.modpow(challenge, &data.key.n()));
            let _1 = &private_commitment.alpha + (challenge * &pdata.plaintext);
            let _3 = &private_commitment.gamma + (challenge * &private_commitment.mu);
            Proof { _1, _2, _3 }
        }

        fn verify(aux: &Aux, data: &Data, commitment: &Commitment, challenge: &Hash, proof: &Proof) -> Result<(), &'static str> {
            let check1 = || {
                let pt = &proof._1 % data.key.n();
                match data.key.encrypt(pt.to_bytes(), Some(proof._2.clone())) {
                    Some((cipher, _nonce)) =>
                        if cipher != commitment.a.modmul(&data.ciphertext.modpow(challenge, data.key.nn()), data.key.nn()) {
                            Err("check1 failed")
                        }
                        else {
                            Ok(())
                        }
                    None => Err("encrypt failed"),
                }
            };
            fn fail_if(b: bool, msg: &'static str) -> Result<(), &'static str>{
                if b {
                    Ok(())
                } else {
                    Err(msg)
                }
            }
            let check2 = |()| fail_if(
                combine(&aux.s, &proof._1, &aux.t, &proof._3, &aux.rsa_modulo)
                == combine(&commitment.c, &1.into(), &commitment.s, challenge, &aux.rsa_modulo),
                "check2",
            );
            let check3 = |()| fail_if(proof._1 <= (BigNumber::one() << (L + EPSILON)), "check3");
            check1().and_then(check2).and_then(check3)
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

    fn combine(l: &BigNumber, le: &BigNumber, r: &BigNumber, re: &BigNumber, m: &BigNumber) -> BigNumber {
        l.modpow(&le, &m).modmul( &r.modpow(&re, &m), &m )
    }

    #[cfg(test)]
    mod test {
        use unknown_order::BigNumber;
        use super::{L, EPSILON};

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

            let ((commitment, challenge, proof), _pcomm) = crate::fiat_shamir::run_scheme::<super::P, _>(&aux, &data, &pdata, rand_core::OsRng::default());
            let r = crate::fiat_shamir::scheme_verify::<super::P>(&aux, &data, &commitment, &challenge, &proof);
            match r {
                Ok(()) => (),
                Err(e) => panic!("{}", e),
            }
        }
        #[test]
        fn failing() {
            let p = BigNumber::from_slice(glass_pumpkin::prime::new(L + 1).unwrap().to_bytes_le());
            let q = BigNumber::from_slice(glass_pumpkin::prime::new(EPSILON + 1).unwrap().to_bytes_le());
            let private_key = libpaillier::DecryptionKey::with_primes_unchecked(&p, &q).unwrap();
            let key = libpaillier::EncryptionKey::from(&private_key);
            let plaintext: BigNumber = BigNumber::one() << (L + EPSILON + 1);
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

            let ((commitment, challenge, proof), _pcomm) = crate::fiat_shamir::run_scheme::<super::P, _>(&aux, &data, &pdata, rand_core::OsRng::default());
            let r = crate::fiat_shamir::scheme_verify::<super::P>(&aux, &data, &commitment, &challenge, &proof);
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
