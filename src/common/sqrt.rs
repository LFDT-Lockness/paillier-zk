use rand_core::RngCore;
use rug::{Complete, Integer};

use crate::IntegerExt;

/// Find principal square root in a Blum modulus quotient ring.
///
/// Pre-requisites:
/// - x is a quadratic residue in Zn
/// - `n = pq`, p and q are Blum primes
/// If these don't hold, the result is a bogus number in Zn
pub fn blum_sqrt(x: &Integer, p: &Integer, q: &Integer, n: &Integer) -> Integer {
    // Exponent in pq Blum modulus to obtain the principal square root.
    // Described in [Handbook of Applied cryptography, p. 75, Fact
    // 2.160](https://cacr.uwaterloo.ca/hac/about/chap2.pdf)
    let e = ((p - 1u8).complete() * (q - 1u8).complete() + 4) / 8;

    // e guaranteed to be non-negative by the prerequisite that p and q are blum primes
    #[allow(clippy::expect_used)]
    x.pow_mod_ref(&e, n)
        .expect("e guaranteed to be non-negative")
        .into()
}

/// Find `(y' = (-1)^a w^b y, a, b)` such that y' is a quadratic residue in Zn.
///
/// a and b are treated as false = 0, true = 1
///
/// Pre-requisites:
/// - `n = pq`, p and q are Blum primes
/// - `jacobi(w, n) = -1`, that is w is quadratic non-residue in Zn with jacobi
/// symbol of -1
/// If these don't hold, the y' might not exist. In this case, returns `None`
pub fn find_residue(
    y: &Integer,
    w: &Integer,
    p: &Integer,
    q: &Integer,
    n: &Integer,
) -> Option<(bool, bool, Integer)> {
    let jp = y.modulo(p).jacobi(p);
    let jq = y.modulo(q).jacobi(q);
    match (jp, jq) {
        (1, 1) => return Some((false, false, y.clone())),
        (-1, -1) => return Some((true, false, (n - y).complete())),
        _ => (),
    }

    let y = (y * w).complete().modulo(n);
    let jp = y.modulo(p).jacobi(p);
    let jq = y.modulo(q).jacobi(q);
    match (jp, jq) {
        (1, 1) => Some((false, true, y)),
        (-1, -1) => Some((true, true, n - y)),
        _ => None,
    }
}

/// Finds a element in Zn that has jacobi symbol of -1
pub fn sample_neg_jacobi<R: RngCore>(n: &Integer, rng: &mut R) -> Integer {
    loop {
        let w = n
            .clone()
            .random_below(&mut fast_paillier::utils::external_rand(rng));
        if w.jacobi(n) == -1 {
            break w;
        }
    }
}
