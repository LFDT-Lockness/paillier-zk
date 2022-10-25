use rand_core::RngCore;
use unknown_order::BigNumber;

/// n = pq
pub fn blum_sqrt(x: &BigNumber, p: &BigNumber, q: &BigNumber, n: &BigNumber) -> BigNumber {
    let e = ((p - 1) * (q - 1) + 4) / 8;
    x.modpow(&e, &n)
}

pub fn find_residue(y: &BigNumber, w: &BigNumber, p: &BigNumber, q: &BigNumber, n: &BigNumber) -> (bool, bool, BigNumber) {
    eprintln!("y = {}\nw = {}\np = {}\nq = {}\nn = {}", y, w, p, q, n);
    for (a, b) in TWO_BOOLS {
        let y = if b { w * y } else { y.clone() };
        let y = if a { n - y } else { y };
        let jp = jacobi(&y, p);
        let jq = jacobi(&y, q);
        if jp == 1 && jq == 1 {
            return (a, b, y)
        } else {
            eprintln!("y = {}\njp = {}, jq = {}", y, jp, jq);
        }
    }
    panic!("so w should have had jacobi of -1, not just be a non-residue")
}

pub const TWO_BOOLS: [(bool, bool); 4] = [
    (false, false),
    (true, false),
    (false, true),
    (true, true),
];


pub fn non_residue_in<R: RngCore>(n: &BigNumber, mut rng: R) -> BigNumber {
    loop {
        let w = BigNumber::from_rng(&n, &mut rng);
        if jacobi(&w, &n) == -1 {
            break w;
        }
    }
}

#[allow(clippy::many_single_char_names)]
pub fn jacobi(x: &BigNumber, y: &BigNumber) -> isize {
    let five = BigNumber::from(5);
    let four = BigNumber::from(4);
    let three = BigNumber::from(3);
    let two = BigNumber::from(2);
    let one = BigNumber::one();
    let zero = BigNumber::zero();
    if y % &two == zero || y <= &zero {
        panic!(
            "invalid arguments, y must be an odd integer,but got {:?}",
            y
        );
    }

    let mut k = y.clone();
    let mut n = x % &k;
    let mut j = 1;

    while n != zero {
        // reduce two's power
        while &n % &two == zero {
            n = &n >> 1;
            let r = &k % BigNumber::from(8);
            if r == three || r == five {
                j = -j;
            }
        }
        core::mem::swap(&mut n, &mut k);
        if &n % &four == three && &k % &four == three {
            j = -j;
        }
        n = n % &k;
    }

    if k == one {
        j
    } else {
        0
    }
}
