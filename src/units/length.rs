use num_traits::{Float, FromPrimitive};

use crate::{label::UnitLabel, private::conversion_constant};

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

impl LengthUnit {
    pub fn scale_to_centimetres<T>(self) -> T
    where
        T: Float + FromPrimitive,
    {
        match self {
            Self::Metre => conversion_constant(100.0),
            Self::Centimetre => T::one(),
            Self::Millimetre => conversion_constant(0.1),
            Self::Micrometre => conversion_constant(1.0e-4),
            Self::Nanometre => conversion_constant(1.0e-7),
        }
    }
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

impl LengthUnit {
    pub fn to_centimetres_factor<T>(self) -> T
    where
        T: Float + FromPrimitive,
    {
        self.scale_to_centimetres()
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
        assert_close(LengthUnit::Metre.scale_to_centimetres::<f64>(), 100.0);
        assert_close(LengthUnit::Centimetre.scale_to_centimetres::<f64>(), 1.0);
        assert_close(LengthUnit::Millimetre.scale_to_centimetres::<f64>(), 1.0e-1);
        assert_close(LengthUnit::Micrometre.scale_to_centimetres::<f64>(), 1.0e-4);
        assert_close(LengthUnit::Nanometre.scale_to_centimetres::<f64>(), 1.0e-7);
    }
}
