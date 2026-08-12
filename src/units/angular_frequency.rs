use num_traits::{Float, FromPrimitive};

use crate::{UnitLabel, private::conversion_constant};

/// Units of angular frequency.
///
/// Angular frequency, denoted by `ω`, is measured in radians per second and is
/// related to ordinary frequency by
///
/// ```text
/// ω = 2πf.
/// ```
///
/// These units provide linear scaling into radians per second. They do not
/// independently identify a canonical TMM coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AngularFrequencyUnit {
    /// Radians per second.
    RadiansPerSecond,

    /// Kiloradians per second.
    KiloradiansPerSecond,

    /// Megaradians per second.
    MegaradiansPerSecond,

    /// Gigaradians per second.
    GigaradiansPerSecond,

    /// Teraradians per second.
    TeraradiansPerSecond,

    /// Petaradians per second.
    PetaradiansPerSecond,
}

impl AngularFrequencyUnit {
    pub fn scale_to_radians_per_second<T>(self) -> T
    where
        T: Float + FromPrimitive,
    {
        match self {
            Self::RadiansPerSecond => T::one(),
            Self::KiloradiansPerSecond => conversion_constant(1e3),
            Self::MegaradiansPerSecond => conversion_constant(1e6),
            Self::GigaradiansPerSecond => conversion_constant(1e9),
            Self::TeraradiansPerSecond => conversion_constant(1e12),
            Self::PetaradiansPerSecond => conversion_constant(1e15),
        }
    }
}

impl UnitLabel for AngularFrequencyUnit {
    fn symbol(&self) -> &'static str {
        match self {
            Self::RadiansPerSecond => "rad·s⁻¹",
            Self::KiloradiansPerSecond => "krad·s⁻¹",
            Self::MegaradiansPerSecond => "Mrad·s⁻¹",
            Self::GigaradiansPerSecond => "Grad·s⁻¹",
            Self::TeraradiansPerSecond => "Trad·s⁻¹",
            Self::PetaradiansPerSecond => "Prad·s⁻¹",
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

    fn angular_frequency_scale(unit: AngularFrequencyUnit) -> f64 {
        unit.scale_to_radians_per_second::<f64>()
    }

    #[test]
    fn angular_frequency_scales_are_correct() {
        let cases = [
            (AngularFrequencyUnit::RadiansPerSecond, 1.0),
            (AngularFrequencyUnit::KiloradiansPerSecond, 1.0e3),
            (AngularFrequencyUnit::MegaradiansPerSecond, 1.0e6),
            (AngularFrequencyUnit::GigaradiansPerSecond, 1.0e9),
            (AngularFrequencyUnit::TeraradiansPerSecond, 1.0e12),
            (AngularFrequencyUnit::PetaradiansPerSecond, 1.0e15),
        ];

        for (unit, expected) in cases {
            assert_close(angular_frequency_scale(unit), expected);
        }
    }

    #[test]
    fn angular_frequency_unit_symbols_are_correct() {
        assert_eq!(AngularFrequencyUnit::RadiansPerSecond.symbol(), "rad·s⁻¹",);
        assert_eq!(
            AngularFrequencyUnit::KiloradiansPerSecond.symbol(),
            "krad·s⁻¹",
        );
        assert_eq!(
            AngularFrequencyUnit::MegaradiansPerSecond.symbol(),
            "Mrad·s⁻¹",
        );
        assert_eq!(
            AngularFrequencyUnit::GigaradiansPerSecond.symbol(),
            "Grad·s⁻¹",
        );
        assert_eq!(
            AngularFrequencyUnit::TeraradiansPerSecond.symbol(),
            "Trad·s⁻¹",
        );
        assert_eq!(
            AngularFrequencyUnit::PetaradiansPerSecond.symbol(),
            "Prad·s⁻¹",
        );
    }
}
