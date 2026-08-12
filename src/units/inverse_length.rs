use crate::{label::UnitLabel, private::conversion_constant};

use num_traits::{Float, FromPrimitive};

/// Units for quantities with dimensions of inverse length.
///
/// Numeric values are scaled internally into inverse centimetres.
///
/// ```text
/// q_cm⁻¹ = q_local × scale
/// ```
///
/// For example:
///
/// ```text
/// 1 m⁻¹  = 10⁻² cm⁻¹
/// 1 µm⁻¹ = 10⁴ cm⁻¹
/// ```
///
/// This type describes only the unit of an inverse-length quantity. It does
/// not identify the physical coordinate represented by that quantity. In
/// particular, it may be used as part of a vacuum-wavenumber, parallel-
/// wavenumber, or other inverse-length coordinate definition.
///
/// [`UnitScale`]: crate::transform::UnitScale
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum InverseLengthUnit {
    /// Inverse metres.
    PerMetre,

    /// Inverse centimetres.
    PerCentimetre,

    /// Inverse millimetres.
    PerMillimetre,

    /// Inverse micrometres.
    PerMicrometre,

    /// Inverse nanometres.
    PerNanometre,
}

impl InverseLengthUnit {
    pub fn scale_to_inverse_centimetres<T>(self) -> T
    where
        T: Float + FromPrimitive,
    {
        match self {
            Self::PerMetre => conversion_constant(1.0e-2),
            Self::PerCentimetre => T::one(),
            Self::PerMillimetre => conversion_constant(10.0),
            Self::PerMicrometre => conversion_constant(1.0e4),
            Self::PerNanometre => conversion_constant(1.0e7),
        }
    }
}

impl UnitLabel for InverseLengthUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::PerMetre => "m⁻¹",
            Self::PerCentimetre => "cm⁻¹",
            Self::PerMillimetre => "mm⁻¹",
            Self::PerMicrometre => "µm⁻¹",
            Self::PerNanometre => "nm⁻¹",
        }
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
    fn inverse_metre_scale_is_correct() {
        let scale = InverseLengthUnit::PerMetre.scale_to_inverse_centimetres::<f64>();

        assert_close(scale, 1.0e-2);
    }

    #[test]
    fn inverse_centimetre_scale_is_correct() {
        let scale = InverseLengthUnit::PerCentimetre.scale_to_inverse_centimetres::<f64>();

        assert_close(scale, 1.0);
    }

    #[test]
    fn inverse_millimetre_scale_is_correct() {
        let scale = InverseLengthUnit::PerMillimetre.scale_to_inverse_centimetres::<f64>();

        assert_close(scale, 10.0);
    }

    #[test]
    fn inverse_micrometre_scale_is_correct() {
        let scale = InverseLengthUnit::PerMicrometre.scale_to_inverse_centimetres::<f64>();

        assert_close(scale, 1.0e4);
    }

    #[test]
    fn inverse_nanometre_scale_is_correct() {
        let scale = InverseLengthUnit::PerNanometre.scale_to_inverse_centimetres::<f64>();

        assert_close(scale, 1.0e7);
    }

    #[test]
    fn scales_convert_values_to_inverse_centimetres() {
        let cases = [
            (InverseLengthUnit::PerMetre, 500.0, 5.0),
            (InverseLengthUnit::PerCentimetre, 500.0, 500.0),
            (InverseLengthUnit::PerMillimetre, 500.0, 5_000.0),
            (InverseLengthUnit::PerMicrometre, 0.5, 5_000.0),
            (InverseLengthUnit::PerNanometre, 5.0e-4, 5_000.0),
        ];

        for (unit, value, expected) in cases {
            let scale = unit.scale_to_inverse_centimetres::<f64>();

            assert_close(value * scale, expected);
        }
    }

    #[test]
    fn scales_convert_inverse_centimetres_back_to_local_units() {
        let canonical = 5_000.0_f64;

        let cases = [
            (InverseLengthUnit::PerMetre, 5.0e5),
            (InverseLengthUnit::PerCentimetre, 5_000.0),
            (InverseLengthUnit::PerMillimetre, 500.0),
            (InverseLengthUnit::PerMicrometre, 0.5),
            (InverseLengthUnit::PerNanometre, 5.0e-4),
        ];

        for (unit, expected) in cases {
            let scale = unit.scale_to_inverse_centimetres::<f64>();

            assert_close(canonical / scale, expected);
        }
    }

    #[test]
    fn all_scales_are_positive_and_finite() {
        let units = [
            InverseLengthUnit::PerMetre,
            InverseLengthUnit::PerCentimetre,
            InverseLengthUnit::PerMillimetre,
            InverseLengthUnit::PerMicrometre,
            InverseLengthUnit::PerNanometre,
        ];

        for unit in units {
            let scale = unit.scale_to_inverse_centimetres::<f64>();

            assert!(scale.is_finite());
            assert!(scale > 0.0);
        }
    }

    #[test]
    fn wavenumber_unit_symbols_are_correct() {
        assert_eq!(InverseLengthUnit::PerMetre.symbol(), "m⁻¹");
        assert_eq!(InverseLengthUnit::PerCentimetre.symbol(), "cm⁻¹");
        assert_eq!(InverseLengthUnit::PerMillimetre.symbol(), "mm⁻¹");
        assert_eq!(InverseLengthUnit::PerMicrometre.symbol(), "µm⁻¹");
        assert_eq!(InverseLengthUnit::PerNanometre.symbol(), "nm⁻¹");
    }

    #[test]
    fn scales_support_f32() {
        let micrometre_scale =
            InverseLengthUnit::PerMicrometre.scale_to_inverse_centimetres::<f32>();

        let nanometre_scale = InverseLengthUnit::PerNanometre.scale_to_inverse_centimetres::<f32>();

        assert!((micrometre_scale - 1.0e4).abs() < 1.0e-2);
        assert!((nanometre_scale - 1.0e7).abs() < 1.0);
    }
}
