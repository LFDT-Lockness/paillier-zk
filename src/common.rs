pub mod sqrt;

use unknown_order::BigNumber;

/// Generate element in Zm*. Does so by trial.
pub fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, mut rng: R) -> BigNumber {
    loop {
        let r = BigNumber::from_rng(modulo, &mut rng);
        if r.gcd(modulo) == 1.into() {
            break r;
        }
    }
}

/// Compute l^le * r^re modulo m
pub fn combine(
    l: &BigNumber,
    le: &BigNumber,
    r: &BigNumber,
    re: &BigNumber,
    m: &BigNumber,
) -> BigNumber {
    l.modpow(&le, &m).modmul(&r.modpow(&re, &m), &m)
}
