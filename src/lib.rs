pub mod paillier_encryption_in_range;
pub mod paillier_affine_operation_in_range;
pub mod paillier_blum_modulus;
mod sqrt;

use unknown_order::BigNumber;

pub mod fiat_shamir {
    use rand_core::RngCore;

    pub trait FiatShamir {
        type Data;
        type PrivateData;
        type Commitment;
        type PrivateCommitment;
        type Hash;
        type Proof;
        type Aux;

        fn commit<R: RngCore>(
            aux: &Self::Aux,
            data: &Self::Data,
            private_data: &Self::PrivateData,
            rng: R,
        ) -> (Self::Commitment, Self::PrivateCommitment);
        fn challenge(
            aux: &Self::Aux,
            data: &Self::Data,
            commitment: &Self::Commitment,
        ) -> Self::Hash;
        fn prove(
            aux: &Self::Aux,
            data: &Self::Data,
            private_data: &Self::PrivateData,
            commitment: &Self::Commitment,
            private_commitment: &Self::PrivateCommitment,
            challenge: &Self::Hash,
        ) -> Self::Proof;
        fn verify(
            aux: &Self::Aux,
            data: &Self::Data,
            commitment: &Self::Commitment,
            challenge: &Self::Hash,
            proof: &Self::Proof,
        ) -> Result<(), &'static str>;
    }

    pub fn run_scheme<F: FiatShamir, R: RngCore>(
        aux: &F::Aux,
        data: &F::Data,
        private_data: &F::PrivateData,
        mut rng: R,
    ) -> ((F::Commitment, F::Hash, F::Proof), F::PrivateCommitment) {
        let (commitment, private_commitment) = F::commit(aux, data, private_data, &mut rng);
        let challenge = F::challenge(aux, data, &commitment);
        let proof = F::prove(
            aux,
            data,
            private_data,
            &commitment,
            &private_commitment,
            &challenge,
        );
        ((commitment, challenge, proof), private_commitment)
    }

    pub fn scheme_verify<F>(
        aux: &F::Aux,
        data: &F::Data,
        commitment: &F::Commitment,
        challenge: &F::Hash,
        proof: &F::Proof,
    ) -> Result<(), &'static str>
    where
        F: FiatShamir,
        F::Hash: PartialEq,
    {
        let challenge2 = F::challenge(aux, data, commitment);
        if *challenge != challenge2 {
            Err("challenge doesn't match")
        } else {
            F::verify(aux, data, commitment, challenge, proof)
        }
    }
}

const L: usize = 228;
const EPSILON: usize = 322;


fn gen_inversible<R: rand_core::RngCore>(modulo: &BigNumber, mut rng: R) -> BigNumber {
    for _ in 0..100 {
        let r = BigNumber::from_rng(modulo, &mut rng);
        if r.gcd(modulo) == 1.into() {
            return r;
        }
    }
    panic!("Attempts exceeded when generating inversible number");
}


fn combine(
    l: &BigNumber,
    le: &BigNumber,
    r: &BigNumber,
    re: &BigNumber,
    m: &BigNumber,
) -> BigNumber {
    l.modpow(&le, &m).modmul(&r.modpow(&re, &m), &m)
}

fn traverse<'a, I, E>(i: I) -> Result<(), E>
where
    I: IntoIterator<Item = &'a dyn Fn() -> Result<(), E>>,
    E: 'a,
{
    for f in i {
        f()?
    }
    Ok(())
}
