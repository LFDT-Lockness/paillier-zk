use subtle::{Choice, ConstantTimeEq};

#[derive(PartialEq, Eq, Clone, Copy, Hash, PartialOrd, Ord, Default)]
pub struct MillionRing(u64);

impl From<u64> for MillionRing {
    fn from(x: u64) -> Self {
        MillionRing(x)
    }
}

const MODULO: u64 = 1_000_000_007;

impl generic_ec_core::Additive for MillionRing {
    fn add(a: &Self, b: &Self) -> Self {
        ((a.0 + b.0) % MODULO).into()
    }

    fn sub(a: &Self, b: &Self) -> Self {
        ((a.0 + MODULO - b.0) % MODULO).into()
    }

    fn negate(a: &Self) -> Self {
        if a.0 == 0 {
            0u64.into()
        } else {
            (MODULO - a.0).into()
        }
    }
}

impl From<generic_ec_core::CurveGenerator> for MillionRing {
    fn from(_: generic_ec_core::CurveGenerator) -> Self {
        2u64.into()
    }
}

impl generic_ec_core::Zero for MillionRing {
    fn zero() -> Self {
        0u64.into()
    }

    fn is_zero(x: &Self) -> Choice {
        if x.0 == 0 {
            1.into()
        } else {
            0.into()
        }
    }
}

impl zeroize::Zeroize for MillionRing {
    fn zeroize(&mut self) {
        self.0 = 0u64
    }
}

impl generic_ec_core::OnCurve for MillionRing {
    #[inline]
    fn is_on_curve(&self) -> Choice {
        Choice::from(1)
    }
}

impl generic_ec_core::SmallFactor for MillionRing {
    #[inline]
    fn is_torsion_free(&self) -> Choice {
        Choice::from(1)
    }
}

impl ConstantTimeEq for MillionRing {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl subtle::ConditionallySelectable for MillionRing {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        u64::conditional_select(&a.0, &b.0, choice).into()
    }
}

impl generic_ec_core::CompressedEncoding for MillionRing {
    type Bytes = [u8; 8];

    fn to_bytes_compressed(&self) -> Self::Bytes {
        self.0.to_le_bytes()
    }
}

impl generic_ec_core::UncompressedEncoding for MillionRing {
    type Bytes = [u8; 8];

    fn to_bytes_uncompressed(&self) -> Self::Bytes {
        self.0.to_le_bytes()
    }
}

impl generic_ec_core::Decode for MillionRing {
    fn decode(bytes: &[u8]) -> Option<Self> {
        let x = u64::from_be_bytes(bytes.try_into().ok()?);
        Some(MillionRing(x % MODULO))
    }
}

#[derive(PartialEq, Eq, Clone, Copy, Hash, PartialOrd, Ord, Default)]
pub struct Scalar(u64);

impl From<u64> for Scalar {
    fn from(x: u64) -> Self {
        Scalar(x % MODULO)
    }
}

impl From<Scalar> for u64 {
    fn from(x: Scalar) -> Self {
        x.0
    }
}

impl generic_ec_core::Additive for Scalar {
    fn add(a: &Self, b: &Self) -> Self {
        ((a.0 + b.0) % MODULO).into()
    }

    fn sub(a: &Self, b: &Self) -> Self {
        ((a.0 + MODULO - b.0) % MODULO).into()
    }

    fn negate(a: &Self) -> Self {
        if a.0 == 0 {
            0u64.into()
        } else {
            (MODULO - a.0).into()
        }
    }
}

impl generic_ec_core::Multiplicative<Scalar> for Scalar {
    type Output = Self;

    fn mul(a: &Self, b: &Self) -> Self {
        (a.0 * b.0 % MODULO).into()
    }
}

impl generic_ec_core::Multiplicative<MillionRing> for Scalar {
    type Output = MillionRing;

    fn mul(a: &Self, b: &MillionRing) -> MillionRing {
        (a.0 * b.0 % MODULO).into()
    }
}

impl generic_ec_core::Multiplicative<generic_ec_core::CurveGenerator> for Scalar {
    type Output = MillionRing;

    fn mul(a: &Self, _: &generic_ec_core::CurveGenerator) -> MillionRing {
        generic_ec_core::Multiplicative::mul(a, &MillionRing::from(generic_ec_core::CurveGenerator))
    }
}

impl generic_ec_core::Invertible for Scalar {
    fn invert(x: &Self) -> subtle::CtOption<Self> {
        subtle::CtOption::new(*x, x.ct_eq(&1u64.into()))
    }
}

impl generic_ec_core::Zero for Scalar {
    fn zero() -> Self {
        0u64.into()
    }

    fn is_zero(x: &Self) -> Choice {
        x.0.ct_eq(&0)
    }
}

impl generic_ec_core::One for Scalar {
    fn one() -> Self {
        1u64.into()
    }

    fn is_one(x: &Self) -> Choice {
        x.0.ct_eq(&1)
    }
}

impl generic_ec_core::Samplable for Scalar {
    fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        Scalar(rng.next_u64())
    }
}

impl zeroize::Zeroize for Scalar {
    fn zeroize(&mut self) {
        self.0 = 0
    }
}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl subtle::ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        u64::conditional_select(&a.0, &b.0, choice).into()
    }
}

impl generic_ec_core::IntegerEncoding for Scalar {
    type Bytes = [u8; 8];

    fn to_be_bytes(&self) -> Self::Bytes {
        self.0.to_be_bytes()
    }

    fn to_le_bytes(&self) -> Self::Bytes {
        self.0.to_le_bytes()
    }

    fn from_be_bytes(bytes: &Self::Bytes) -> Self {
        u64::from_be_bytes(*bytes).into()
    }

    fn from_le_bytes(bytes: &Self::Bytes) -> Self {
        u64::from_le_bytes(*bytes).into()
    }

    fn from_be_bytes_exact(bytes: &Self::Bytes) -> Option<Self> {
        Some(u64::from_be_bytes(*bytes).into())
    }

    fn from_le_bytes_exact(bytes: &Self::Bytes) -> Option<Self> {
        Some(u64::from_le_bytes(*bytes).into())
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, PartialOrd, Ord, Default)]
pub struct C;

impl generic_ec::Curve for C {
    const CURVE_NAME: &'static str = "test curve Z/1000000007";
    type Point = MillionRing;
    type Scalar = Scalar;
    type CompressedPointArray = [u8; 8];
    type UncompressedPointArray = [u8; 8];
    type ScalarArray = [u8; 8];
    type CoordinateArray = [u8; 8];
}

impl generic_ec_core::hash_to_curve::HashToCurve for C {
    fn hash_to_curve(
        tag: generic_ec::hash_to_curve::Tag,
        msgs: &[&[u8]],
    ) -> Result<Self::Point, generic_ec_core::Error> {
        use sha2::Digest;
        let mut digest = sha2::Sha256::new();
        digest.update(tag.as_bytes());
        for msg in msgs {
            digest.update(msg);
        }
        let bytes = digest.finalize();
        let bytes = bytes.as_slice()[0..8].try_into().unwrap();
        let x = u64::from_be_bytes(bytes) % MODULO;
        Ok(MillionRing(x))
    }

    fn hash_to_scalar(
        tag: generic_ec::hash_to_curve::Tag,
        msgs: &[&[u8]],
    ) -> Result<Self::Scalar, generic_ec_core::Error> {
        use sha2::Digest;
        let mut digest = sha2::Sha256::new();
        digest.update(tag.as_bytes());
        for msg in msgs {
            digest.update(msg);
        }
        let bytes = digest.finalize();
        let bytes = bytes.as_slice()[0..8].try_into().unwrap();
        let x = u64::from_be_bytes(bytes);
        Ok(Scalar(x))
    }
}
