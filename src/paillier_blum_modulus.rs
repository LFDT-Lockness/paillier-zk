use rand_core::RngCore;
use unknown_order::BigNumber;

use crate::sqrt::*;

pub struct P;

const M: usize = 13;

pub struct Data {
    n: BigNumber,
}

pub struct PrivateData {
    p: BigNumber,
    q: BigNumber,
}

pub struct Commitment {
    w: BigNumber,
}

pub type PrivateCommitment = ();

#[derive(Debug, PartialEq, Eq)]
pub struct Hash {
    ys: Vec<BigNumber>,
}

pub struct ProofPoint {
    x: BigNumber,
    a: bool,
    b: bool,
    z: BigNumber,
}

pub struct Proof {
    points: Vec<ProofPoint>,
}

pub type Aux = ();

impl crate::fiat_shamir::FiatShamir for P {
    type Data = Data;
    type PrivateData = PrivateData;
    type Commitment = Commitment;
    type PrivateCommitment = PrivateCommitment;
    type Hash = Hash;
    type Proof = Proof;
    type Aux = Aux;

    fn commit<R: RngCore>(
        _aux: &Aux,
        Data { ref n }: &Data,
        _pdata: &PrivateData,
        rng: R,
    ) -> (Commitment, ()) {
        /*
        let w = loop {
            // FIXME so which jacobi should it be?
            // Judging by notes in taurusgroup, it should be non-residue in Zn
            // with any jacobi symbol
            let r = BigNumber::from_rng(&n, &mut rng);
            if jacobi(&r, &pdata.p) == -1 && jacobi(&r, &pdata.q) == -1 {
                break r;
            }
        };
        */
        let w = non_residue_in(&n, rng);
        let r = Commitment {
            w
        };
        (r, ())
    }

    fn challenge(_aux: &Aux, _data: &Data, _commitment: &Commitment) -> Hash {
        let mut ys = Vec::new();
        ys.reserve_exact(M);
        for i in 1..=M {
            ys.push(i.into());
        }
        Hash { ys }
    }

    fn prove(
        (): &Aux,
        Data { ref n }: &Data,
        PrivateData { ref p, ref q }: &PrivateData,
        Commitment { ref w }: &Commitment,
        (): &PrivateCommitment,
        challenge: &Hash,
    ) -> Self::Proof {
        let sqrt = |x| {
            blum_sqrt(&x, &p, &q, &n)
        };
        let phi = (p - 1) * (q - 1);
        let n_inverse = n.extended_gcd(&phi).x;
        assert_eq!(n_inverse.modmul(&(n % &phi), &phi), BigNumber::one());
        let points = challenge
            .ys
            .iter()
            .map(|y| {
                let z = y.modpow(&n_inverse, n);
                let (a, b, y_) = find_residue(&y, w, p, q, n);
                let x = sqrt(sqrt(y_));
                ProofPoint { x, a, b, z }
            })
            .collect();
        Proof { points }
    }

    fn verify(
        _aux: &Aux,
        data: &Data,
        commitment: &Commitment,
        challenge: &Hash,
        proof: &Proof,
    ) -> Result<(), &'static str> {
        if data.n.is_prime() {
            return Err("N is prime");
        }
        if &data.n % BigNumber::from(2) == BigNumber::zero() {
            return Err("N is even");
        }
        if proof.points.len() != challenge.ys.len() {
            return Err("Incorrect amount of proofs");
        }
        for (point, y) in proof.points.iter().zip(challenge.ys.iter()) {
            if point.z.modpow(&data.n, &data.n) != *y {
                return Err("z^n != y");
            }
            let m1 = if point.a {
                &data.n - BigNumber::one()
            } else {
                BigNumber::one()
            };
            let one = BigNumber::one();
            let m2 = if point.b { &commitment.w } else { &one };
            let rhs = m1.modmul(m2, &data.n).modmul(y, &data.n);
            if point.x.modpow(&4.into(), &data.n) != rhs {
                return Err("x^4 != rhs");
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use unknown_order::BigNumber;

    #[test]
    fn passing() {
        let mut rng = rand_core::OsRng::default();
        let p = blum_prime(256);
        let q = blum_prime(256);
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let ((commitment, challenge, proof), _pcomm) =
            crate::fiat_shamir::run_scheme::<super::P, _>(
                &(),
                &data,
                &pdata,
                &mut rng,
            );
        let r = crate::fiat_shamir::scheme_verify::<super::P>(
            &(),
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
        let mut rng = rand_core::OsRng::default();
        let p = BigNumber::prime(256);
        let q = loop { // non blum prime
            let q = BigNumber::prime(256);
            if &q % 4 == BigNumber::one() {
                break q;
            }
        };
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let ((commitment, challenge, proof), _pcomm) =
            crate::fiat_shamir::run_scheme::<super::P, _>(
                &(),
                &data,
                &pdata,
                &mut rng,
            );
        let r = crate::fiat_shamir::scheme_verify::<super::P>(
            &(),
            &data,
            &commitment,
            &challenge,
            &proof,
        );
        match r {
            Ok(()) => panic!("should have failed"),
            Err(_) => (),
        }
    }

    fn blum_prime(s: usize) -> BigNumber {
        let three = BigNumber::from(3);
        loop {
            let p = BigNumber::prime(s);
            if &p % 4 == three {
                break p;
            }
        }
    }
}
