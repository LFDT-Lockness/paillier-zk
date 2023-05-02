use sha2::digest;
use sha2::digest::Digest;

pub struct HashRng<F, D: Digest> {
    hash: F,
    counter: u64,
    buffer: digest::Output<D>,
    offset: usize,
}

impl<F, D: Digest> HashRng<F, D>
{
    pub fn new(hash: F) -> Self
    where
        F: Fn(D) -> digest::Output<D>,
    {
        let d: D = D::new().chain_update(0u64.to_le_bytes());
        let buffer: digest::Output<D> = hash(d);
        HashRng {
            hash,
            counter: 1,
            offset: 0,
            buffer,
        }
    }
}

impl<F, D> rand_core::RngCore for HashRng<F, D>
where
    D: Digest,
    digest::Output<D>: std::borrow::Borrow<[u8]>,
    F: Fn(D) -> digest::Output<D>,
{
    fn next_u32(&mut self) -> u32 {
        use std::borrow::Borrow;

        const SIZE: usize = std::mem::size_of::<u32>();
        if self.offset + SIZE > self.buffer.borrow().len() {
            self.buffer = (self.hash) (D::new().chain_update(self.counter.to_le_bytes()));
            self.counter += 1;
            self.offset = 0;
        }
        let bytes = &self.buffer.borrow()[self.offset .. self.offset + SIZE];
        self.offset += SIZE;
        let bytes: [u8; SIZE] = bytes.try_into().unwrap();
        u32::from_le_bytes(bytes)
    }

    fn next_u64(&mut self) -> u64 {
        rand_core::impls::next_u64_via_u32(self)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        rand_core::impls::fill_bytes_via_next(self, dest)
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        Ok(self.fill_bytes(dest))
    }
}

#[cfg(test)]
mod test {
    use rand_core::RngCore;
    use sha2::Digest;

    #[test]
    fn generate_bytes() {
        let hash = |d: sha2::Sha256| d
            .chain_update("foobar")
            .finalize();
        let mut rng = super::HashRng::new(hash);

        let mut zeroes = 0;
        let mut total = 0;

        let mut bytes = [0; 32];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 128];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 1];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 1];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 32];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 128];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);

        let mut bytes = [0; 137];
        rng.fill_bytes(&mut bytes);
        total += bytes.len();
        zeroes += bytes.iter().filter(|b| **b == 0).count();
        assert!(zeroes <= total + 256 / 256);
    }
}
