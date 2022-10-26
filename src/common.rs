use unknown_order::BigNumber;

/// Generate element in Zm*. Does so by trial.
pub fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, mut rng: R) -> BigNumber {
    for _ in 0..100 {
        let r = BigNumber::from_rng(modulo, &mut rng);
        if r.gcd(modulo) == 1.into() {
            return r;
        }
    }
    panic!("Attempts exceeded when generating inversible number");
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

/// Collect results from iterator, stopping at first error
pub fn sequence<'a, I, E>(i: I) -> Result<(), E>
where
    I: IntoIterator<Item = &'a dyn Fn() -> Result<(), E>>,
    E: 'a,
{
    for f in i {
        f()?
    }
    Ok(())
}
