//! ZK-proof of Paillier-Blum modulus. Called Пmod or Rmod in the CGGMP21 paper.
//!
//! ## Description
//! A party P has a modulus `N = pq`, with p and q being Blum primes, and
//! `gcd(N, phi(N)) = 1`. P wants to prove that those equalities about N hold,
//! without disclosing p and q.
//!
//! ## Example
//! 0. Prover P derives two Blum primes and makes a Paillier-Blum modulus
//!     ``` no_run
//!     # use paillier_zk::unknown_order::BigNumber;
//!     fn blum_prime(s: usize) -> BigNumber {
//!         let three = BigNumber::from(3);
//!         loop {
//!             let p = BigNumber::prime(s);
//!             if &p % 4 == three {
//!                 break p;
//!             }
//!         }
//!     }
//!     let p = blum_prime(256);
//!     let q = blum_prime(256);
//!     let n = &p * &q;
//!     // Prover can then make a key from it
//!     let pkey = libpaillier::DecryptionKey::with_primes_unchecked(&p, &q);
//!     ```
//! 1. P computes a non-interactive proof that `n` is a Paillier-Blum modulus:
//!     ``` no_run
//!     use paillier_zk::paillier_blum_modulus as p;
//!     # use generic_ec_core::hash_to_curve::Tag;
//!     # let (n, p, q) = todo!();
//!     const TAG: Tag = Tag::new_unwrap("application name".as_bytes());
//!     const SECURITY: usize = 33;
//!
//!     let data = p::Data { n };
//!     let pdata = p::PrivateData { p, q };
//!     let mut rng = rand_core::OsRng::default();
//!
//!     let (commitment, challenge, proof) =
//!         p::compute_proof::<{SECURITY}, _>(
//!             TAG,
//!             &data,
//!             &pdata,
//!             &mut rng,
//!         );
//!     ```
//! 2. P sends `data, commitment, challenge, proof` to the verifier V
//! 3. V verifies the proof:
//!     ``` no_run
//!     # use paillier_zk::paillier_blum_modulus as p;
//!     # let (data, commitment, challenge, proof) = todo!();
//!     # const SECURITY: usize = 33;
//!     p::verify::<{SECURITY}>(
//!         &data,
//!         &commitment,
//!         &challenge,
//!         &proof,
//!     );
//!     ```
//! 4. If the verification succeeded, V can continue communication with P

use crate::unknown_order::BigNumber;
use generic_ec::hash_to_curve::Tag;
use rand_core::RngCore;

use crate::common::sqrt::{blum_sqrt, find_residue, non_residue_in};

#[derive(Debug, PartialEq, Eq)]
pub enum InvalidProof {
    ModulusIsPrime,
    ModulusIsEven,
    IncorrectNthRoot,
    IncorrectFourthRoot,
}

/// Public data that both parties know: the Paillier-Blum modulus
pub struct Data {
    pub n: BigNumber,
}

/// Private data of prover
pub struct PrivateData {
    pub p: BigNumber,
    pub q: BigNumber,
}

/// Prover's first message, obtained by `commit`
pub struct Commitment {
    pub w: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// `challenge`
///
/// Consists of `M` singular challenges
#[derive(Debug, PartialEq, Eq)]
pub struct Challenge<const M: usize> {
    pub ys: [BigNumber; M],
}

pub struct ProofPoint {
    pub x: BigNumber,
    pub a: bool,
    pub b: bool,
    pub z: BigNumber,
}

/// The ZK proof. Computed by `prove`. Consists of M proofs for each challenge
pub struct Proof<const M: usize> {
    pub points: [ProofPoint; M],
}

/// Create random commitment
pub fn commit<R: RngCore>(Data { ref n }: &Data, rng: R) -> Commitment {
    Commitment {
        w: non_residue_in(n, rng),
    }
}

/// Deterministically compute challenge based on prior known values in protocol
pub fn challenge<const M: usize>(tag: Tag, data: &Data, commitment: &Commitment) -> Challenge<M> {
    // since we can't use Default and BigNumber isn't copy, we initialize
    // like this
    let mut ys = [(); M].map(|()| BigNumber::zero());
    for (i, y_ref) in ys.iter_mut().enumerate() {
        *y_ref = crate::common::hash2field::hash_to_field(
            tag,
            &data.n,
            &[
                &data.n.to_bytes(),
                &commitment.w.to_bytes(),
                &(i as u64 + 1).to_le_bytes(),
            ],
        );
    }
    Challenge { ys }
}

/// Compute proof for given data and prior protocol values
pub fn prove<const M: usize>(
    Data { ref n }: &Data,
    PrivateData { ref p, ref q }: &PrivateData,
    Commitment { ref w }: &Commitment,
    challenge: &Challenge<M>,
) -> Proof<M> {
    let sqrt = |x| blum_sqrt(&x, p, q, n);
    let phi = (p - 1) * (q - 1);
    let n_inverse = n.extended_gcd(&phi).x;
    assert_eq!(n_inverse.modmul(&(n % &phi), &phi), BigNumber::one());
    let points = challenge.ys.clone().map(|y| {
        let z = y.modpow(&n_inverse, n);
        let (a, b, y_) = find_residue(&y, w, p, q, n);
        let x = sqrt(sqrt(y_));
        ProofPoint { x, a, b, z }
    });
    Proof { points }
}

/// Verify the proof. If this succeeds, the relation Rmod holds with chance
/// `1/2^M`
pub fn verify<const M: usize>(
    data: &Data,
    commitment: &Commitment,
    challenge: &Challenge<M>,
    proof: &Proof<M>,
) -> Result<(), InvalidProof> {
    if data.n.is_prime() {
        return Err(InvalidProof::ModulusIsPrime);
    }
    if &data.n % BigNumber::from(2) == BigNumber::zero() {
        return Err(InvalidProof::ModulusIsEven);
    }
    for (point, y) in proof.points.iter().zip(challenge.ys.iter()) {
        if point.z.modpow(&data.n, &data.n) != *y {
            return Err(InvalidProof::IncorrectNthRoot);
        }
        let y = y.clone();
        let y = if point.a { &data.n - y } else { y };
        let y = if point.b {
            y.modmul(&commitment.w, &data.n)
        } else {
            y
        };
        if point.x.modpow(&4.into(), &data.n) != y {
            return Err(InvalidProof::IncorrectFourthRoot);
        }
    }
    Ok(())
}

/// Compute proof for the given data, producing random commitment and
/// deriving determenistic challenge.
///
/// Obtained from the above interactive proof via Fiat-Shamir heuristic.
pub fn compute_proof<const M: usize, R: RngCore>(
    tag: Tag,
    data: &Data,
    pdata: &PrivateData,
    rng: R,
) -> (Commitment, Challenge<M>, Proof<M>) {
    let commitment = commit(data, rng);
    let challenge = challenge(tag, data, &commitment);
    let proof = prove(data, pdata, &commitment, &challenge);
    (commitment, challenge, proof)
}

#[cfg(test)]
mod test {
    use crate::unknown_order::BigNumber;

    #[test]
    fn passing() {
        let mut rng = rand_core::OsRng::default();
        let p = blum_prime(256);
        let q = blum_prime(256);
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) =
            super::compute_proof::<65, _>(tag, &data, &pdata, &mut rng);
        let r = super::verify(&data, &commitment, &challenge, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }

    #[test]
    fn failing() {
        let mut rng = rand_core::OsRng::default();
        let p = BigNumber::prime(256);
        let q = loop {
            // non blum prime
            let q = BigNumber::prime(256);
            if &q % 4 == BigNumber::one() {
                break q;
            }
        };
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let (commitment, challenge, proof) =
            super::compute_proof::<65, _>(tag, &data, &pdata, &mut rng);
        let r = super::verify(&data, &commitment, &challenge, &proof);
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
