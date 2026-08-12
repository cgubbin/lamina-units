use crate::{label::UnitLabel, private::conversion_constant};

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

impl AngleUnit {
    pub fn scale_to_radians<T: Float + FromPrimitive>(&self) -> T {
        match self {
            Self::Radian => T::one(),
            Self::Degree => conversion_constant(std::f64::consts::PI / 180.0),
        }
    }
}

impl UnitLabel for AngleUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::Radian => "rad",
            Self::Degree => "°",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::label::UnitLabel;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-14_f64.max(expected.abs() * 1.0e-13);

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}",
        );
    }

    #[test]
    fn radian_scale_is_identity() {
        let scale = AngleUnit::Radian.scale_to_radians::<f64>();

        assert_eq!(scale, 1.0);
    }

    #[test]
    fn degree_scale_converts_to_radians() {
        let scale = AngleUnit::Degree.scale_to_radians::<f64>();

        assert_close(scale, std::f64::consts::PI / 180.0);
    }

    #[test]
    fn one_hundred_and_eighty_degrees_is_pi_radians() {
        let scale = AngleUnit::Degree.scale_to_radians::<f64>();

        assert_close(180.0 * scale, std::f64::consts::PI);
    }

    #[test]
    fn ninety_degrees_is_pi_over_two_radians() {
        let scale = AngleUnit::Degree.scale_to_radians::<f64>();

        assert_close(90.0 * scale, std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn angle_symbols_are_correct() {
        assert_eq!(AngleUnit::Radian.symbol(), "rad");
        assert_eq!(AngleUnit::Degree.symbol(), "°");
    }

    #[test]
    fn angle_scale_supports_f32() {
        let scale = AngleUnit::Degree.scale_to_radians::<f32>();

        let expected = std::f32::consts::PI / 180.0;

        assert!((scale - expected).abs() < 1.0e-7);
    }
}
