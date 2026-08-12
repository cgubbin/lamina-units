use std::fmt;

use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{
    AngularFrequencyUnit, Frequency, FrequencyUnit, UnitLabel, VacuumAngularWavenumber, Wavelength,
    private::two_pi,
};

/// Angular frequency `ω`.
///
/// Angular frequency is related to ordinary frequency by
///
/// ```text
/// ω = 2π f.
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AngularFrequency<F> {
    value: F,
    unit: AngularFrequencyUnit,
}

impl<F> AngularFrequency<F> {
    /// Construct an angular frequency.
    pub const fn new(value: F, unit: AngularFrequencyUnit) -> Self {
        Self { value, unit }
    }

    pub fn radians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::RadiansPerSecond)
    }

    pub fn kiloradians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::KiloradiansPerSecond)
    }

    pub fn megaradians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::MegaradiansPerSecond)
    }

    pub fn gigaradians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::GigaradiansPerSecond)
    }

    pub fn teraradians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::TeraradiansPerSecond)
    }

    pub fn petaradians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::PetaradiansPerSecond)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> AngularFrequencyUnit {
        self.unit
    }

    pub fn into_parts(self) -> (F, AngularFrequencyUnit) {
        (self.value, self.unit)
    }

    /// Return the angular frequency expressed in radians per second.
    pub fn as_radians_per_second(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_radians_per_second::<F>() * self.value
    }

    /// Convert to ordinary frequency.
    pub fn frequency(&self) -> Frequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        Frequency::new(
            self.as_radians_per_second() / two_pi::<F>(),
            FrequencyUnit::Hertz,
        )
    }

    /// Convert to vacuum angular wavenumber.
    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        self.frequency().vacuum_angular_wavenumber()
    }

    /// Convert to wavelength
    pub fn wavelength(&self) -> Wavelength<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        self.frequency().wavelength()
    }
}

impl<F> fmt::Display for AngularFrequency<F>
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
        let frequency = AngularFrequency::teraradians_per_second(2.5);

        assert_eq!(frequency.value(), 2.5);
        assert_eq!(frequency.unit(), AngularFrequencyUnit::TeraradiansPerSecond,);
    }

    #[test]
    fn equivalent_units_have_equal_radians_per_second_values() {
        let values = [
            AngularFrequency::radians_per_second(1.0e12),
            AngularFrequency::kiloradians_per_second(1.0e9),
            AngularFrequency::megaradians_per_second(1.0e6),
            AngularFrequency::gigaradians_per_second(1.0e3),
            AngularFrequency::teraradians_per_second(1.0),
        ];

        for value in values {
            assert_eq!(value.as_radians_per_second(), 1.0e12,);
        }
    }

    #[test]
    fn angular_frequency_converts_to_frequency() {
        let omega = AngularFrequency::radians_per_second(std::f64::consts::TAU);

        assert_eq!(omega.frequency().as_hertz(), 1.0,);
    }

    #[test]
    fn frequency_round_trip_is_consistent() {
        let original = Frequency::terahertz(30.0);

        let round_trip = original.angular_frequency().frequency();

        assert_eq!(round_trip.as_hertz(), original.as_hertz(),);
    }

    #[test]
    fn angular_frequency_converts_to_vacuum_angular_wavenumber() {
        let omega = AngularFrequency::radians_per_second(
            std::f64::consts::TAU * crate::constants::SPEED_OF_LIGHT_CM_PER_SECOND,
        );

        assert_eq!(
            omega.vacuum_angular_wavenumber().as_inverse_centimetres(),
            std::f64::consts::TAU,
        );
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(
            AngularFrequency::teraradians_per_second(5.0).to_string(),
            "5 Trad·s⁻¹",
        );
    }
}
