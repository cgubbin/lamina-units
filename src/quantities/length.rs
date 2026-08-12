use num_traits::{Float, FromPrimitive, Zero};
use std::fmt;

use crate::{LengthUnit, UnitLabel};

/// A physical length with an explicit unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Length<F> {
    value: F,
    unit: LengthUnit,
}

impl<F> Length<F> {
    /// Construct a physical length.
    pub const fn new(value: F, unit: LengthUnit) -> Self {
        Self { value, unit }
    }

    /// Construct a length in centimetres.
    pub fn centimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Centimetre)
    }

    /// Construct a length in millimetres.
    pub fn millimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Millimetre)
    }

    /// Construct a length in micrometres.
    pub fn micrometres(value: F) -> Self {
        Self::new(value, LengthUnit::Micrometre)
    }

    /// Construct a length in nanometres.
    pub fn nanometres(value: F) -> Self {
        Self::new(value, LengthUnit::Nanometre)
    }

    /// Return the numerical value in the stored unit.
    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    /// Return the stored unit.
    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    /// Construct a zero length.
    pub fn zero() -> Self
    where
        F: Zero,
    {
        Self {
            value: F::zero(),
            unit: LengthUnit::Centimetre,
        }
    }

    /// Return half this length, preserving its unit.
    pub fn half(self) -> Self
    where
        F: Float,
    {
        self.scale_by(F::one() / (F::one() + F::one()))
    }

    /// Scale this length by a dimensionless factor, preserving its unit.
    pub fn scale_by(self, scalar: F) -> Self
    where
        F: Float,
    {
        let (value, unit) = self.into_parts();

        Self::new(value * scalar, unit)
    }

    /// Return whether the stored value is zero.
    pub fn is_zero(&self) -> bool
    where
        F: Float,
    {
        self.value == F::zero()
    }

    /// Consume this length and return its value and unit.
    pub fn into_parts(self) -> (F, LengthUnit) {
        (self.value, self.unit)
    }

    /// Return the length expressed in centimetres.
    pub fn as_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_centimetres::<F>() * self.value
    }
}

impl<F> fmt::Display for Length<F>
where
    F: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_preserve_value_and_unit() {
        let cases = [
            (Length::centimetres(2.0), 2.0, LengthUnit::Centimetre),
            (Length::millimetres(2.0), 2.0, LengthUnit::Millimetre),
            (Length::micrometres(2.0), 2.0, LengthUnit::Micrometre),
            (Length::nanometres(2.0), 2.0, LengthUnit::Nanometre),
        ];

        for (length, value, unit) in cases {
            assert_eq!(length.value(), value);
            assert_eq!(length.unit(), unit);
        }
    }

    #[test]
    fn equivalent_lengths_have_equal_centimetre_values() {
        let lengths = [
            Length::centimetres(1.0),
            Length::millimetres(10.0),
            Length::micrometres(10_000.0),
            Length::nanometres(10_000_000.0),
        ];

        for length in lengths {
            assert_eq!(length.as_centimetres(), 1.0);
        }
    }

    #[test]
    fn scaling_preserves_unit() {
        let length = Length::nanometres(500.0).scale_by(2.0);

        assert_eq!(length.value(), 1000.0);
        assert_eq!(length.unit(), LengthUnit::Nanometre);
    }

    #[test]
    fn half_preserves_unit() {
        let length = Length::micrometres(2.0).half();

        assert_eq!(length.value(), 1.0);
        assert_eq!(length.unit(), LengthUnit::Micrometre);
    }

    #[test]
    fn zero_is_zero() {
        let length = Length::<f64>::zero();

        assert!(length.is_zero());
        assert_eq!(length.as_centimetres(), 0.0);
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(Length::nanometres(500.0).to_string(), "500 nm",);
    }
}
