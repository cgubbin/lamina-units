use num_traits::FromPrimitive;

pub(crate) fn conversion_constant<T>(value: f64) -> T
where
    T: FromPrimitive,
{
    T::from_f64(value).expect("unit conversion constant should be representable")
}
