use std::fmt;

use num_traits::{Float, FromPrimitive};

use crate::{AngleUnit, UnitLabel};

/// A physical angle with an explicit angular unit.
///
/// `Angle` represents an angle independently of any optical interpretation.
/// In particular, it does not imply an incidence angle or perform any
/// conversion to an in-plane optical coordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Angle<F> {
    value: F,
    unit: AngleUnit,
}

impl<F> Angle<F> {
    /// Construct an angle.
    pub const fn new(value: F, unit: AngleUnit) -> Self {
        Self { value, unit }
    }

    /// Construct an angle in radians.
    pub fn radians(value: F) -> Self {
        Self::new(value, AngleUnit::Radian)
    }

    /// Construct an angle in degrees.
    pub fn degrees(value: F) -> Self {
        Self::new(value, AngleUnit::Degree)
    }

    /// Return the numerical value in the stored unit.
    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    /// Return the stored angular unit.
    pub const fn unit(&self) -> AngleUnit {
        self.unit
    }

    /// Consume this angle and return its value and unit.
    pub fn into_parts(self) -> (F, AngleUnit) {
        (self.value, self.unit)
    }

    /// Return the angle expressed in radians.
    pub fn as_radians(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_radians::<F>() * self.value
    }
}

impl<F> fmt::Display for Angle<F>
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
        let degrees = Angle::degrees(45.0);

        assert_eq!(degrees.value(), 45.0);
        assert_eq!(degrees.unit(), AngleUnit::Degree);

        let radians = Angle::radians(1.5);

        assert_eq!(radians.value(), 1.5);
        assert_eq!(radians.unit(), AngleUnit::Radian);
    }

    #[test]
    fn degrees_convert_to_radians() {
        let angle = Angle::degrees(180.0);

        assert_eq!(angle.as_radians(), std::f64::consts::PI,);
    }

    #[test]
    fn equivalent_angles_have_equal_radian_values() {
        let degrees = Angle::degrees(45.0);
        let radians = Angle::radians(std::f64::consts::FRAC_PI_4);

        assert_eq!(degrees.as_radians(), radians.as_radians(),);
    }

    #[test]
    fn into_parts_preserves_value_and_unit() {
        let angle = Angle::degrees(30.0);

        assert_eq!(angle.into_parts(), (30.0, AngleUnit::Degree),);
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(Angle::degrees(45.0).to_string(), "45 °",);
    }
}
