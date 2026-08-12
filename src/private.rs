use num_traits::{Float, FloatConst, FromPrimitive};

pub(crate) fn conversion_constant<T>(value: f64) -> T
where
    T: FromPrimitive,
{
    T::from_f64(value).expect("unit conversion constant should be representable")
}

pub(crate) fn two_pi<F>() -> F
where
    F: FloatConst + Float,
{
    F::PI() + F::PI()
}
