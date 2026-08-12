use crate::{label::UnitLabel, private::conversion_constant};

use num_traits::{Float, FromPrimitive};

/// Units of ordinary frequency.
///
/// Ordinary frequency, conventionally denoted by `f` or `ν`, counts cycles
/// per unit time. It is related to angular frequency by
///
/// ```text
/// ω = 2πf.
/// ```
///
/// These units provide linear scaling into hertz. They do not themselves
/// define a canonical TMM coordinate; that conversion is performed by a
/// spectral-coordinate type or by [`frequency_to_vacuum_angular_wavenumber`].
///
/// [`frequency_to_vacuum_angular_wavenumber`]:
///     crate::frequency::frequency_to_vacuum_angular_wavenumber
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FrequencyUnit {
    /// Hertz.
    Hertz,

    /// Kilohertz.
    Kilohertz,

    /// Megahertz.
    Megahertz,

    /// Gigahertz.
    Gigahertz,

    /// Terahertz.
    Terahertz,

    /// Petahertz.
    Petahertz,
}

impl FrequencyUnit {
    pub fn scale_to_hertz<T>(self) -> T
    where
        T: Float + FromPrimitive,
    {
        match self {
            Self::Hertz => T::one(),
            Self::Kilohertz => conversion_constant(1.0e3),
            Self::Megahertz => conversion_constant(1.0e6),
            Self::Gigahertz => conversion_constant(1.0e9),
            Self::Terahertz => conversion_constant(1.0e12),
            Self::Petahertz => conversion_constant(1.0e15),
        }
    }
}

impl UnitLabel for FrequencyUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::Hertz => "Hz",
            Self::Kilohertz => "kHz",
            Self::Megahertz => "MHz",
            Self::Gigahertz => "GHz",
            Self::Terahertz => "THz",
            Self::Petahertz => "PHz",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-14_f64.max(1.0e-13 * expected.abs());

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}; \
             absolute error was {:e}",
            (actual - expected).abs(),
        );
    }

    fn frequency_scale(unit: FrequencyUnit) -> f64 {
        unit.scale_to_hertz::<f64>()
    }

    #[test]
    fn ordinary_frequency_scales_are_correct() {
        let cases = [
            (FrequencyUnit::Hertz, 1.0),
            (FrequencyUnit::Kilohertz, 1.0e3),
            (FrequencyUnit::Megahertz, 1.0e6),
            (FrequencyUnit::Gigahertz, 1.0e9),
            (FrequencyUnit::Terahertz, 1.0e12),
            (FrequencyUnit::Petahertz, 1.0e15),
        ];

        for (unit, expected) in cases {
            assert_close(frequency_scale(unit), expected);
        }
    }

    #[test]
    fn ordinary_frequency_unit_symbols_are_correct() {
        assert_eq!(FrequencyUnit::Hertz.symbol(), "Hz");
        assert_eq!(FrequencyUnit::Kilohertz.symbol(), "kHz");
        assert_eq!(FrequencyUnit::Megahertz.symbol(), "MHz");
        assert_eq!(FrequencyUnit::Gigahertz.symbol(), "GHz");
        assert_eq!(FrequencyUnit::Terahertz.symbol(), "THz");
        assert_eq!(FrequencyUnit::Petahertz.symbol(), "PHz");
    }
}
