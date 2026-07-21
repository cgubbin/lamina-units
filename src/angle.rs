use crate::{label::UnitLabel, private::conversion_constant, transform::UnitScale};

use num_traits::{Float, FromPrimitive};

/// Units used to express an angle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AngleUnit {
    /// Radians.
    Radian,

    /// Degrees.
    Degree,
}

impl UnitLabel for AngleUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::Radian => "rad",
            Self::Degree => "°",
        }
    }
}

impl<T> UnitScale<T> for AngleUnit
where
    T: Float + FromPrimitive,
{
    /// Returns the scale from the represented angle unit to radians.
    fn scale(self) -> T {
        match self {
            Self::Radian => T::one(),
            Self::Degree => conversion_constant(std::f64::consts::PI / 180.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::{label::UnitLabel, transform::UnitScale};

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-14_f64.max(expected.abs() * 1.0e-13);

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}",
        );
    }

    #[test]
    fn radian_scale_is_identity() {
        let scale = <AngleUnit as UnitScale<f64>>::scale(AngleUnit::Radian);

        assert_eq!(scale, 1.0);
    }

    #[test]
    fn degree_scale_converts_to_radians() {
        let scale = <AngleUnit as UnitScale<f64>>::scale(AngleUnit::Degree);

        assert_close(scale, std::f64::consts::PI / 180.0);
    }

    #[test]
    fn one_hundred_and_eighty_degrees_is_pi_radians() {
        let scale = <AngleUnit as UnitScale<f64>>::scale(AngleUnit::Degree);

        assert_close(180.0 * scale, std::f64::consts::PI);
    }

    #[test]
    fn ninety_degrees_is_pi_over_two_radians() {
        let scale = <AngleUnit as UnitScale<f64>>::scale(AngleUnit::Degree);

        assert_close(90.0 * scale, std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn angle_symbols_are_correct() {
        assert_eq!(AngleUnit::Radian.symbol(), "rad");
        assert_eq!(AngleUnit::Degree.symbol(), "°");
    }

    #[test]
    fn angle_scale_supports_f32() {
        let scale = <AngleUnit as UnitScale<f32>>::scale(AngleUnit::Degree);

        let expected = std::f32::consts::PI / 180.0;

        assert!((scale - expected).abs() < 1.0e-7);
    }
}
