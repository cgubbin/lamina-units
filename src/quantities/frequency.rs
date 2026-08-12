use std::fmt;

use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{
    AngularFrequency, FrequencyUnit, UnitLabel, VacuumAngularWavenumber, VacuumWavenumber,
    Wavelength, constants::SPEED_OF_LIGHT_CM_PER_SECOND,
};

/// Ordinary frequency `f`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Frequency<F> {
    value: F,
    unit: FrequencyUnit,
}

impl<F> Frequency<F> {
    pub const fn new(value: F, unit: FrequencyUnit) -> Self {
        Self { value, unit }
    }

    pub fn hertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Hertz)
    }

    pub fn kilohertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Kilohertz)
    }

    pub fn megahertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Megahertz)
    }

    pub fn gigahertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Gigahertz)
    }

    pub fn terahertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Terahertz)
    }

    pub fn petahertz(value: F) -> Self {
        Self::new(value, FrequencyUnit::Petahertz)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> FrequencyUnit {
        self.unit
    }

    pub fn into_parts(self) -> (F, FrequencyUnit) {
        (self.value, self.unit)
    }

    /// Return the frequency expressed in hertz.
    pub fn as_hertz(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_hertz::<F>() * self.value
    }

    /// Convert to spectroscopic vacuum wavenumber.
    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FromPrimitive,
    {
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND)
            .expect("speed of light should be representable");

        VacuumWavenumber::inverse_centimetres(self.as_hertz() / c)
    }

    /// Convert to vacuum angular wavenumber
    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        self.vacuum_wavenumber().vacuum_angular_wavenumber()
    }

    /// Convert to angular frequency.
    pub fn angular_frequency(&self) -> AngularFrequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        AngularFrequency::radians_per_second(crate::private::two_pi::<F>() * self.as_hertz())
    }

    /// Convert to wavelength
    pub fn wavelength(&self) -> Wavelength<F>
    where
        F: Float + FromPrimitive,
    {
        self.vacuum_wavenumber().wavelength()
    }
}

impl<F> fmt::Display for Frequency<F>
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
        let frequency = Frequency::terahertz(30.0);

        assert_eq!(frequency.value(), 30.0);
        assert_eq!(frequency.unit(), FrequencyUnit::Terahertz);
    }

    #[test]
    fn equivalent_units_have_equal_hertz_values() {
        let values = [
            Frequency::hertz(1.0e12),
            Frequency::kilohertz(1.0e9),
            Frequency::megahertz(1.0e6),
            Frequency::gigahertz(1.0e3),
            Frequency::terahertz(1.0),
        ];

        for value in values {
            assert_eq!(value.as_hertz(), 1.0e12);
        }
    }

    #[test]
    fn frequency_converts_to_vacuum_wavenumber() {
        let frequency = Frequency::hertz(SPEED_OF_LIGHT_CM_PER_SECOND);

        assert_eq!(frequency.vacuum_wavenumber().as_inverse_centimetres(), 1.0,);
    }

    #[test]
    fn frequency_converts_to_vacuum_angular_wavenumber() {
        let frequency = Frequency::hertz(SPEED_OF_LIGHT_CM_PER_SECOND);

        let k0 = frequency.vacuum_angular_wavenumber();

        assert_eq!(k0.as_inverse_centimetres(), std::f64::consts::TAU,);
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(Frequency::terahertz(30.0).to_string(), "30 THz",);
    }
}
