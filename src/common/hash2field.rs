//! Hashing to field as described in
//! https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11
use generic_ec::hash_to_curve::Tag;
use sha2::digest::generic_array::GenericArray;

use crate::unknown_order::BigNumber;

// Maximum possible security value
const SECURITY_FOR_SHA256: usize = 16;

/// Hash the given messages and domain tag, and produce a number modulo `q`. The
/// bytes generated are guaranteed to be uniform, and the number is guaranteed
/// to be statistically close to uniform if `q` is prime
pub fn hash_to_field(dst: Tag, q: &BigNumber, messages: &[&[u8]]) -> BigNumber {
    // count = 1, m = 1 (since q is prime)
    // L = ceil((ceil(log2(p)) + k) / 8) where k is security parameter, and then
    // len_in_bytes = count * m * L = L
    // In RustCrypto k = 0, but here I chose the maximum possible value; for
    // BigNumber unlike crypto-bigint it's still fast enough and without size
    // problems.
    let len_in_bytes = q.to_bytes().len() + SECURITY_FOR_SHA256;
    let bytes = expand_message_xmd(messages, len_in_bytes, dst);
    BigNumber::from_slice(bytes) % q
}

/// The expand_message_xmd function produces a uniformly random byte string
/// using a cryptographic hash function H that outputs b bits.
/// This function should follow some requirements, see here:
/// https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-5.4.1
fn expand_message_xmd(messages: &[&[u8]], len_in_bytes: usize, dst: Tag) -> Vec<u8> {
    use sha2::Digest;
    // It seems that security properties of this function are independent of
    // choice of k, at least if 2k <= b - bits produced by hash function.

    let b_in_bytes = sha2::Sha256::output_size();
    // in the rfc, it's a byte, but I think it's unnecessary for our case
    let blocks = (len_in_bytes + b_in_bytes - 1) / b_in_bytes;

    let mut b_0 = sha2::Sha256::new();
    // update with msg_prime
    // skip padding
    for message in messages {
        b_0.update(message);
    }
    b_0.update((len_in_bytes as u16).to_le_bytes());
    b_0.update([0]);
    b_0.update(dst.as_bytes());
    b_0.update((dst.as_bytes().len() as u64).to_le_bytes());
    let b_0 = b_0.finalize();

    let mut result = Vec::new();
    result.reserve(blocks * b_in_bytes);

    let mut b_cur = sha2::Sha256::new();
    b_cur.update(b_0);
    b_cur.update([1]);
    b_cur.update(dst.as_bytes());
    b_cur.update((dst.as_bytes().len() as u64).to_le_bytes());
    let mut b_prev = b_cur.finalize();
    result.extend_from_slice(&b_prev);
    for i in 2..=blocks {
        // this assumes that block size is less than 2^64, i.e. len_in_bytes <
        // 2^72, which I think is fair
        let i = i as u64;
        let mut b_cur = sha2::Sha256::new();
        b_cur.update(strxor(b_prev, b_0));
        b_cur.update(i.to_le_bytes());
        b_cur.update(dst.as_bytes());
        b_cur.update((dst.as_bytes().len() as u64).to_le_bytes());
        b_prev = b_cur.finalize();
        result.extend_from_slice(&b_prev);
    }

    result
}

fn strxor<L: sha2::digest::generic_array::ArrayLength<u8>>(
    arr1: GenericArray<u8, L>,
    arr2: GenericArray<u8, L>,
) -> GenericArray<u8, L> {
    let mut r: GenericArray<u8, L> = Default::default();
    for i in 0..L::to_usize() {
        r[i] = arr1[i] ^ arr2[i];
    }
    r
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
        let message = "1";
        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());
        let r = super::hash_to_field(tag, &q, &[message.as_bytes()]);
        assert!(r < q, "result too big");
        let lower_bound = q / 2;
        assert!(r > lower_bound, "result too small");

        // small results very well can happen, but at least for "1" the result
        // is big, which means big results happen as well
    }
}
