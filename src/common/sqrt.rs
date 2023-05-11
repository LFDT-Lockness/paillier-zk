use crate::unknown_order::BigNumber;
use rand_core::RngCore;

/// Find principal square root in a Blum modulus quotient ring.
///
/// Pre-requisites:
/// - x is a quadratic residue in Zn
/// - `n = pq`, p and q are Blum primes
/// If these don't hold, the result is a bogus number in Zn
pub fn blum_sqrt(x: &BigNumber, p: &BigNumber, q: &BigNumber, n: &BigNumber) -> BigNumber {
    // Exponent in pq Blum modulus to obtain the principal square root.
    // Described in [Handbook of Applied cryptography, p. 75, Fact
    // 2.160](https://cacr.uwaterloo.ca/hac/about/chap2.pdf)
    let e = ((p - 1) * (q - 1) + 4) / 8;

    // e guaranteed to be non-negative by the prerequisite that p and q are blum primes
    #[allow(clippy::disallowed_methods)]
    x.modpow(&e, n)
}

/// Find `(y' = (-1)^a w^b y, a, b)` such that y' is a quadratic residue in Zn.
///
/// a and b are treated as false = 0, true = 1
///
/// Pre-requisites:
/// - `n = pq`, p and q are Blum primes
/// - `jacobi(w, n) = -1`, that is w is quadratic non-residue in Zn with jacobi
/// symbol of -1
/// If these don't hold, the y' might not exist. In this case, returns `y' = 0`
pub fn find_residue(
    y: &BigNumber,
    w: &BigNumber,
    p: &BigNumber,
    q: &BigNumber,
    n: &BigNumber,
) -> (bool, bool, BigNumber) {
    for (a, b) in TWO_BOOLS {
        let y = if b { w.modmul(y, n) } else { y.clone() };
        let y = if a { n - y } else { y };
        let jp = jacobi(&(&y % p), p);
        let jq = jacobi(&(&y % q), q);
        if jp == 1 && jq == 1 {
            return (a, b, y);
        }
    }
    (false, false, BigNumber::zero())
}

const TWO_BOOLS: [(bool, bool); 4] = [(false, false), (true, false), (false, true), (true, true)];

/// Find a quadratic non-residue in Zn. Does so by generating a random number
/// and checking its jacobi symbol. The return value is guaranteed to have
/// jacobi symbol of -1
pub fn non_residue_in<R: RngCore>(n: &BigNumber, mut rng: R) -> BigNumber {
    loop {
        let w = BigNumber::from_rng(n, &mut rng);
        if jacobi(&w, n) == -1 {
            break w;
        }
    }
}

/// Compute jacobi symbol of a over n
///
/// Requires odd `n >= 3`, and `0 <= a < n`. If it doesn't hold, function panics if debug asserts are enabled,
/// or returns 0 otherwise.
///
/// Implementation is taken from [Handbook of Applied cryptography][book], p. 75, Algorithm 2.149
///
/// [book]: https://cacr.uwaterloo.ca/hac/about/chap2.pdf
#[inline(always)]
pub fn jacobi(a: &BigNumber, n: &BigNumber) -> isize {
    jacobi_inner(1, a, n)
}

/// Computes jacobi symbol of `a` over `n` multiplied at `mult`
///
/// `mult` is a small modification over original algorithm defined in the book, it helps to keep
/// function in tail recursion form, which ensures that recursion is optimized out
#[allow(clippy::if_same_then_else, clippy::identity_op)]
fn jacobi_inner(mult: isize, a: &BigNumber, n: &BigNumber) -> isize {
    let one = &BigNumber::one();
    let two = &BigNumber::from(2);
    let three = &BigNumber::from(3);

    // Validate inputs
    if !(n % two == *one && n >= three && BigNumber::zero() <= *a && a < n) {
        debug_assert!(false, "invalid inputs: a = {a}, n = {n}");
        return 0;
    }

    // Step 1
    if a.is_zero() {
        return 0;
    }
    // Step 2
    if a.is_one() {
        return mult * 1;
    }

    // Step 3. Find a1, e such that a = 2^e * a1 where a1 is odd
    let mut a1 = a.clone();
    let mut e = 0;
    while &a1 % two != *one {
        e += 1;
        a1 = a1 >> 1;
    }
    debug_assert_eq!(*a, &a1 * (one << e));

    // Step 4
    let n_mod_8 = n % 8;
    let mut s = if e % 2 == 0 {
        // if e is even, s = 1
        1
    } else if n_mod_8 == *one || n_mod_8 == BigNumber::from(7) {
        // if n = 1 or 7 (mod 8), s = 1
        1
    } else if n_mod_8 == *three || n_mod_8 == BigNumber::from(5) {
        // if n = 3 or 5 (mod 8), s = -1
        -1
    } else {
        unreachable!()
    };

    // Step 5
    if n % 4 == *three && &a1 % 4 == *three {
        s = -s
    }

    // Step 6
    let n1 = n % &a1;

    // Step 7
    if a1.is_one() {
        mult * s
    } else {
        jacobi_inner(mult * s, &n1, &a1)
    }
}

#[cfg(test)]
mod test {
    use crate::unknown_order::BigNumber;

    #[test]
    fn jacobi_and_sqrt() {
        let mut rng = rand_dev::DevRng::new();
        // Create a blum modulus for the calculations to make sense
        let p = BigNumber::safe_prime_from_rng(128, &mut rng);
        let q = BigNumber::safe_prime_from_rng(128, &mut rng);
        let n = &p * &q;

        for _ in 0..100 {
            let x = BigNumber::from_rng(&n, &mut rng);
            let root = super::blum_sqrt(&x, &p, &q, &n);
            let x_ = root.modmul(&root, &n);

            let jp = super::jacobi(&(&x % &p), &p);
            let jq = super::jacobi(&(&x % &q), &q);
            let j = super::jacobi(&x, &n);
            assert_eq!(jp * jq, j);

            if jp == 1 && jq == 1 {
                assert_eq!(x_, x);
            } else {
                assert_ne!(x_, x);
            }
        }
    }
}
