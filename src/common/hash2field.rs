//! Hashing to field with a prng
use crate::unknown_order::BigNumber;
use generic_ec::hash_to_curve::Tag;

/// Hash the given messages and domain tag, and produce a number modulo `q`. The
/// bytes generated are guaranteed to be uniform, and the number is guaranteed
/// to be statistically close to uniform if `q` is prime
pub fn hash_to_field(dst: Tag, q: &BigNumber, messages: &[&[u8]]) -> BigNumber {
    use rand_core::SeedableRng;
    use sha2::Digest;

    let mut digest = sha2::Sha256::new();
    digest.update(dst.as_bytes());
    digest.update((messages.len() as u64).to_le_bytes());
    for message in messages {
        digest.update(message);
    }
    let mut rng = rand_chacha::ChaCha20Rng::from_seed(digest.finalize().into());
    BigNumber::from_rng(q, &mut rng)
}

#[cfg(test)]
mod test {
    use libpaillier::unknown_order::BigNumber;

    #[test]
    fn hash_to_field() {
        let q = BigNumber::from(1_000_000_007);
        let message = "a message";
        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let r = super::hash_to_field(tag, &q, &[message.as_bytes()]);
        assert!(r < q, "result too big");
        assert_ne!(r, BigNumber::zero());
        assert!(r > BigNumber::zero(), "negative result");
    }

    #[test]
    fn big_hash() {
        let q = BigNumber::from_slice([
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]);
        let message = "2";
        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let r = super::hash_to_field(tag, &q, &[message.as_bytes()]);
        assert!(r < q, "result too big");
        let lower_bound = q / 2;
        assert!(r > lower_bound, "result too small");

        // small results very well can happen, but at least for "2" the result
        // is big, which means big results happen as well
    }
}
