use crate::{
    label::UnitLabel,
    private::conversion_constant,
    transform::{TransformJet, UnitScale},
};

use num_traits::{Float, FromPrimitive};

/// The speed of light in vacuum, in centimetres per second.
///
/// This is the exact SI value
///
/// ```text
/// c = 299_792_458 m·s⁻¹
///   = 2.997_924_58 × 10¹⁰ cm·s⁻¹.
/// ```
pub const SPEED_OF_LIGHT_CM_PER_S: f64 = 2.997_924_58e10;

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

impl<T> UnitScale<T> for FrequencyUnit
where
    T: Float + FromPrimitive,
{
    fn scale(self) -> T {
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

impl<T> UnitScale<T> for AngularFrequencyUnit
where
    T: Float + FromPrimitive,
{
    fn scale(self) -> T {
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

/// Converts ordinary frequency into canonical vacuum angular wavenumber.
///
/// The canonical coordinate is
///
/// ```text
/// k₀ = 2πf / c,
/// ```
///
/// where `f` is converted into hertz and `c` is expressed in centimetres per
/// second. The resulting `k₀` is therefore expressed in inverse centimetres.
///
/// The returned jet contains
///
/// ```text
/// value  = k₀(f)
/// first  = dk₀/df
/// second = d²k₀/df² = 0,
/// ```
///
/// where the derivatives are with respect to the numeric frequency value in
/// the supplied [`FrequencyUnit`].
///
/// This function uses angular wavenumber `k₀ = 2π/λ`, not the spectroscopic
/// wavenumber `1/λ`.
pub(crate) fn frequency_to_vacuum_angular_wavenumber<T>(
    frequency: T,
    unit: FrequencyUnit,
) -> TransformJet<T>
where
    T: Float + FromPrimitive,
{
    let frequency_scale: T = unit.scale();
    let two_pi: T = conversion_constant(std::f64::consts::TAU);
    let c = conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

    let scale = two_pi * frequency_scale / c;

    TransformJet::new(frequency * scale, scale, T::zero())
}

pub(crate) fn vacuum_angular_wavenumber_to_frequency<T>(
    vacuum_angular_wavenumber: T,
    unit: FrequencyUnit,
) -> TransformJet<T>
where
    T: Float + FromPrimitive,
{
    let frequency_scale: T = unit.scale();
    let two_pi: T = conversion_constant(std::f64::consts::TAU);
    let c: T = conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

    let scale = c / (two_pi * frequency_scale);

    TransformJet::new(vacuum_angular_wavenumber * scale, scale, T::zero())
}

/// Converts angular frequency into canonical vacuum angular wavenumber.
///
/// The canonical coordinate is
///
/// ```text
/// k₀ = ω / c,
/// ```
///
/// where `ω` is converted into radians per second and `c` is expressed in
/// centimetres per second. The resulting `k₀` is expressed in inverse
/// centimetres.
///
/// The returned jet contains
///
/// ```text
/// value  = k₀(ω)
/// first  = dk₀/dω
/// second = d²k₀/dω² = 0,
/// ```
///
/// where the derivatives are with respect to the numeric angular-frequency
/// value in the supplied [`AngularFrequencyUnit`].
pub(crate) fn angular_frequency_to_vacuum_angular_wavenumber<T>(
    angular_frequency: T,
    unit: AngularFrequencyUnit,
) -> TransformJet<T>
where
    T: Float + FromPrimitive,
{
    let angular_frequency_scale: T = unit.scale();
    let scale = angular_frequency_scale / conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

    TransformJet::new(angular_frequency * scale, scale, T::zero())
}

pub(crate) fn vacuum_angular_wavenumber_to_angular_frequency<T>(
    vacuum_angular_wavenumber: T,
    unit: AngularFrequencyUnit,
) -> TransformJet<T>
where
    T: Float + FromPrimitive,
{
    let angular_frequency_scale: T = unit.scale();
    let c: T = conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

    let scale = c / angular_frequency_scale;

    TransformJet::new(vacuum_angular_wavenumber * scale, scale, T::zero())
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
        <FrequencyUnit as UnitScale<f64>>::scale(unit)
    }

    fn angular_frequency_scale(unit: AngularFrequencyUnit) -> f64 {
        <AngularFrequencyUnit as UnitScale<f64>>::scale(unit)
    }

    #[test]
    fn speed_of_light_has_expected_value() {
        assert_eq!(SPEED_OF_LIGHT_CM_PER_S, 2.997_924_58e10);
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
    fn ordinary_frequency_unit_symbols_are_correct() {
        assert_eq!(FrequencyUnit::Hertz.symbol(), "Hz");
        assert_eq!(FrequencyUnit::Kilohertz.symbol(), "kHz");
        assert_eq!(FrequencyUnit::Megahertz.symbol(), "MHz");
        assert_eq!(FrequencyUnit::Gigahertz.symbol(), "GHz");
        assert_eq!(FrequencyUnit::Terahertz.symbol(), "THz");
        assert_eq!(FrequencyUnit::Petahertz.symbol(), "PHz");
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

    #[test]
    fn one_hertz_converts_to_expected_vacuum_angular_wavenumber() {
        let transform = frequency_to_vacuum_angular_wavenumber(1.0_f64, FrequencyUnit::Hertz);

        let expected = std::f64::consts::TAU / SPEED_OF_LIGHT_CM_PER_S;

        assert_close(transform.value, expected);
        assert_close(transform.first, expected);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn one_radian_per_second_converts_to_expected_vacuum_angular_wavenumber() {
        let transform = angular_frequency_to_vacuum_angular_wavenumber(
            1.0_f64,
            AngularFrequencyUnit::RadiansPerSecond,
        );

        let expected = 1.0 / SPEED_OF_LIGHT_CM_PER_S;

        assert_close(transform.value, expected);
        assert_close(transform.first, expected);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn terahertz_converts_to_expected_vacuum_angular_wavenumber() {
        let frequency = 30.0_f64;

        let transform = frequency_to_vacuum_angular_wavenumber(frequency, FrequencyUnit::Terahertz);

        let expected_scale = std::f64::consts::TAU * 1.0e12 / SPEED_OF_LIGHT_CM_PER_S;

        assert_close(transform.value, frequency * expected_scale);
        assert_close(transform.first, expected_scale);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn frequency_and_angular_frequency_are_consistent() {
        let frequency_hz = 30.0e12_f64;
        let angular_frequency = std::f64::consts::TAU * frequency_hz;

        let from_frequency =
            frequency_to_vacuum_angular_wavenumber(frequency_hz, FrequencyUnit::Hertz);

        let from_angular_frequency = angular_frequency_to_vacuum_angular_wavenumber(
            angular_frequency,
            AngularFrequencyUnit::RadiansPerSecond,
        );

        assert_close(from_frequency.value, from_angular_frequency.value);
    }

    #[test]
    fn scaled_frequency_units_are_consistent() {
        let from_hertz = frequency_to_vacuum_angular_wavenumber(30.0e12_f64, FrequencyUnit::Hertz);

        let from_terahertz =
            frequency_to_vacuum_angular_wavenumber(30.0_f64, FrequencyUnit::Terahertz);

        assert_close(from_hertz.value, from_terahertz.value);
    }

    #[test]
    fn scaled_angular_frequency_units_are_consistent() {
        let angular_frequency = 30.0e12_f64;

        let from_radians_per_second = angular_frequency_to_vacuum_angular_wavenumber(
            angular_frequency,
            AngularFrequencyUnit::RadiansPerSecond,
        );

        let from_teraradians_per_second = angular_frequency_to_vacuum_angular_wavenumber(
            30.0_f64,
            AngularFrequencyUnit::TeraradiansPerSecond,
        );

        assert_close(
            from_radians_per_second.value,
            from_teraradians_per_second.value,
        );
    }

    #[test]
    fn frequency_transform_is_linear() {
        let transform = frequency_to_vacuum_angular_wavenumber(12.5_f64, FrequencyUnit::Gigahertz);

        assert_close(transform.value, 12.5 * transform.first);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn angular_frequency_transform_is_linear() {
        let transform = angular_frequency_to_vacuum_angular_wavenumber(
            12.5_f64,
            AngularFrequencyUnit::GigaradiansPerSecond,
        );

        assert_close(transform.value, 12.5 * transform.first);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn frequency_conversion_supports_f32() {
        let transform = frequency_to_vacuum_angular_wavenumber(30.0_f32, FrequencyUnit::Terahertz);

        let expected = std::f32::consts::TAU * 30.0e12 / SPEED_OF_LIGHT_CM_PER_S as f32;

        let relative_error = ((transform.value - expected) / expected).abs();

        assert!(relative_error < 1.0e-6);
        assert_eq!(transform.second, 0.0);
    }

    #[test]
    fn angular_frequency_conversion_supports_f32() {
        let transform = angular_frequency_to_vacuum_angular_wavenumber(
            30.0_f32,
            AngularFrequencyUnit::TeraradiansPerSecond,
        );

        let expected = 30.0e12 / SPEED_OF_LIGHT_CM_PER_S as f32;

        let relative_error = ((transform.value - expected) / expected).abs();

        assert!(relative_error < 1.0e-6);
        assert_eq!(transform.second, 0.0);
    }
}
