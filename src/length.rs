use num_traits::{Float, FromPrimitive};

use crate::{
    label::UnitLabel,
    private::conversion_constant,
    transform::{TransformJet, UnitScale},
};

/// A unit used to express a physical length.
///
/// Lengths are converted to and from the canonical length unit used by the
/// TMM backend, which is the centimetre.
///
/// The associated canonical transformations return a [`TransformJet`] whose
/// derivatives are taken with respect to the numeric value expressed in this
/// unit. For example, for micrometres,
///
/// ```text
/// x_cm = 10⁻⁴ x_µm
/// dx_cm / dx_µm = 10⁻⁴
/// d²x_cm / dx_µm² = 0
/// ```
///
/// [`TransformJet`]: crate::TransformJet
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LengthUnit {
    /// Metres.
    Metre,

    /// Centimetres, the canonical backend length unit.
    Centimetre,

    /// Millimetres.
    Millimetre,

    /// Micrometres.
    Micrometre,

    /// Nanometres.
    Nanometre,
}

impl UnitLabel for LengthUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::Metre => "m",
            Self::Centimetre => "cm",
            Self::Millimetre => "mm",
            Self::Micrometre => "µm",
            Self::Nanometre => "nm",
        }
    }
}

impl<T> UnitScale<T> for LengthUnit
where
    T: Float + FromPrimitive,
{
    fn scale(self) -> T {
        match self {
            Self::Metre => conversion_constant(100.0),
            Self::Centimetre => T::one(),
            Self::Millimetre => conversion_constant(0.1),
            Self::Micrometre => conversion_constant(1.0e-4),
            Self::Nanometre => conversion_constant(1.0e-7),
        }
    }
}

impl LengthUnit {
    pub fn to_canonical<T>(self, value: T) -> TransformJet<T>
    where
        T: Float + FromPrimitive,
    {
        let scale = self.scale();

        TransformJet::new(value * scale, scale, T::zero())
    }

    pub fn from_canonical<T>(self, value: T) -> TransformJet<T>
    where
        T: Float + FromPrimitive,
    {
        let inverse_scale = T::one() / self.scale();

        TransformJet::new(value * inverse_scale, inverse_scale, T::zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-14_f64.max(1.0e-14 * expected.abs());

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}; \
             absolute error was {:e}",
            (actual - expected).abs(),
        );
    }

    #[test]
    fn length_unit_scales_are_correct() {
        assert_close(
            <LengthUnit as UnitScale<f64>>::scale(LengthUnit::Metre),
            100.0,
        );
        assert_close(
            <LengthUnit as UnitScale<f64>>::scale(LengthUnit::Centimetre),
            1.0,
        );
        assert_close(
            <LengthUnit as UnitScale<f64>>::scale(LengthUnit::Millimetre),
            1.0e-1,
        );
        assert_close(
            <LengthUnit as UnitScale<f64>>::scale(LengthUnit::Micrometre),
            1.0e-4,
        );
        assert_close(
            <LengthUnit as UnitScale<f64>>::scale(LengthUnit::Nanometre),
            1.0e-7,
        );
    }

    #[test]
    fn metres_convert_to_centimetres() {
        let transform = LengthUnit::Metre.to_canonical(500.0_f64);

        assert_close(transform.value, 50_000.0);
        assert_close(transform.first, 100.0);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn centimetres_do_not_change_when_converted_to_canonical() {
        let transform = LengthUnit::Centimetre.to_canonical(500.0_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 1.0);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn millimetres_convert_to_centimetres() {
        let transform = LengthUnit::Millimetre.to_canonical(500.0_f64);

        assert_close(transform.value, 50.0);
        assert_close(transform.first, 1.0e-1);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn micrometres_convert_to_centimetres() {
        let transform = LengthUnit::Micrometre.to_canonical(500.0_f64);

        assert_close(transform.value, 5.0e-2);
        assert_close(transform.first, 1.0e-4);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn nanometres_convert_to_centimetres() {
        let transform = LengthUnit::Nanometre.to_canonical(500.0_f64);

        assert_close(transform.value, 5.0e-5);
        assert_close(transform.first, 1.0e-7);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn centimetres_convert_to_metres() {
        let transform = LengthUnit::Metre.from_canonical(50_000.0_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 1.0e-2);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn canonical_centimetres_do_not_change() {
        let transform = LengthUnit::Centimetre.from_canonical(500.0_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 1.0);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn centimetres_convert_to_millimetres() {
        let transform = LengthUnit::Millimetre.from_canonical(50.0_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 10.0);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn centimetres_convert_to_micrometres() {
        let transform = LengthUnit::Micrometre.from_canonical(5.0e-2_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 1.0e4);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn centimetres_convert_to_nanometres() {
        let transform = LengthUnit::Nanometre.from_canonical(5.0e-5_f64);

        assert_close(transform.value, 500.0);
        assert_close(transform.first, 1.0e7);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn all_length_units_round_trip_through_canonical() {
        let units = [
            LengthUnit::Metre,
            LengthUnit::Centimetre,
            LengthUnit::Millimetre,
            LengthUnit::Micrometre,
            LengthUnit::Nanometre,
        ];

        for unit in units {
            let original = 12.345_f64;

            let canonical = unit.to_canonical(original);
            let recovered = unit.from_canonical(canonical.value);

            assert_close(recovered.value, original);
        }
    }

    #[test]
    fn all_length_units_round_trip_from_canonical() {
        let units = [
            LengthUnit::Metre,
            LengthUnit::Centimetre,
            LengthUnit::Millimetre,
            LengthUnit::Micrometre,
            LengthUnit::Nanometre,
        ];

        for unit in units {
            let original_canonical = 12.345_f64;

            let local = unit.from_canonical(original_canonical);
            let recovered = unit.to_canonical(local.value);

            assert_close(recovered.value, original_canonical);
        }
    }

    #[test]
    fn forward_and_inverse_first_derivatives_are_reciprocal() {
        let units = [
            LengthUnit::Metre,
            LengthUnit::Centimetre,
            LengthUnit::Millimetre,
            LengthUnit::Micrometre,
            LengthUnit::Nanometre,
        ];

        for unit in units {
            let forward = unit.to_canonical(12.345_f64);
            let inverse = unit.from_canonical(forward.value);

            assert_close(forward.first * inverse.first, 1.0);
        }
    }

    #[test]
    fn composed_forward_and_inverse_second_derivative_is_zero() {
        let units = [
            LengthUnit::Metre,
            LengthUnit::Centimetre,
            LengthUnit::Millimetre,
            LengthUnit::Micrometre,
            LengthUnit::Nanometre,
        ];

        for unit in units {
            let forward = unit.to_canonical(12.345_f64);
            let inverse = unit.from_canonical(forward.value);

            // For x(q(x)):
            //
            // d²x/dx² =
            //     d²x/dq² (dq/dx)²
            //     + dx/dq d²q/dx².
            let composed_second =
                inverse.second * forward.first * forward.first + inverse.first * forward.second;

            assert_close(composed_second, 0.0);
        }
    }

    #[test]
    fn length_unit_symbols_are_correct() {
        assert_eq!(LengthUnit::Metre.symbol(), "m");
        assert_eq!(LengthUnit::Centimetre.symbol(), "cm");
        assert_eq!(LengthUnit::Millimetre.symbol(), "mm");
        assert_eq!(LengthUnit::Micrometre.symbol(), "µm");
        assert_eq!(LengthUnit::Nanometre.symbol(), "nm");
    }

    #[test]
    fn conversion_supports_f32() {
        let forward = LengthUnit::Micrometre.to_canonical(500.0_f32);
        let inverse = LengthUnit::Micrometre.from_canonical(forward.value);

        assert!((forward.value - 5.0e-2).abs() < 1.0e-6);
        assert!((forward.first - 1.0e-4).abs() < 1.0e-8);
        assert_eq!(forward.second, 0.0);

        assert!((inverse.value - 500.0).abs() < 1.0e-3);
        assert!((inverse.first - 1.0e4).abs() < 1.0);
        assert_eq!(inverse.second, 0.0);
    }
}
