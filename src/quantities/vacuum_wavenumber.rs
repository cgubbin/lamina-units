use std::fmt;

use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{InverseLengthUnit, UnitLabel, VacuumAngularWavenumber, Wavelength};

/// Spectroscopic vacuum wavenumber.
///
/// The vacuum wavenumber is
///
/// ```text
/// ν̃ = 1 / λ₀.
/// ```
///
/// It has dimensions of inverse length and is distinct from the vacuum angular
/// wavenumber `k₀ = 2π / λ₀`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VacuumWavenumber<F> {
    value: F,
    unit: InverseLengthUnit,
}

impl<F> VacuumWavenumber<F> {
    /// Construct a vacuum wavenumber.
    pub const fn new(value: F, unit: InverseLengthUnit) -> Self {
        Self { value, unit }
    }

    pub fn inverse_metres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerMetre)
    }

    pub fn inverse_centimetres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerCentimetre)
    }

    pub fn inverse_millimetres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerMillimetre)
    }

    pub fn inverse_micrometres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerMicrometre)
    }

    pub fn inverse_nanometres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerNanometre)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> InverseLengthUnit {
        self.unit
    }

    pub fn into_parts(self) -> (F, InverseLengthUnit) {
        (self.value, self.unit)
    }

    /// Return the vacuum wavenumber in inverse centimetres.
    pub fn as_inverse_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_inverse_centimetres::<F>() * self.value
    }

    /// Convert to vacuum angular wavenumber.
    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        VacuumAngularWavenumber::inverse_centimetres(
            crate::private::two_pi::<F>() * self.as_inverse_centimetres(),
        )
    }

    /// Convert to vacuum wavelength.
    pub fn wavelength(&self) -> Wavelength<F>
    where
        F: Float + FromPrimitive,
    {
        Wavelength::centimetres(F::one() / self.as_inverse_centimetres())
    }
}

impl<F> fmt::Display for VacuumWavenumber<F>
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
        let wavenumber = VacuumWavenumber::inverse_centimetres(1200.0);

        assert_eq!(wavenumber.value(), 1200.0);
        assert_eq!(wavenumber.unit(), InverseLengthUnit::PerCentimetre,);
    }

    #[test]
    fn equivalent_units_have_equal_inverse_centimetre_values() {
        let values = [
            VacuumWavenumber::inverse_metres(1.0e5),
            VacuumWavenumber::inverse_centimetres(1.0e3),
            VacuumWavenumber::inverse_millimetres(100.0),
            VacuumWavenumber::inverse_micrometres(0.1),
            VacuumWavenumber::inverse_nanometres(1.0e-4),
        ];

        for value in values {
            assert_eq!(value.as_inverse_centimetres(), 1.0e3);
        }
    }

    #[test]
    fn converts_to_vacuum_angular_wavenumber() {
        let wavenumber = VacuumWavenumber::inverse_centimetres(1000.0);

        let k0 = wavenumber.vacuum_angular_wavenumber();

        assert_eq!(
            k0.as_inverse_centimetres(),
            2.0 * std::f64::consts::PI * 1000.0,
        );
    }

    #[test]
    fn round_trip_through_angular_wavenumber() {
        let original = VacuumWavenumber::inverse_centimetres(1234.5);

        let round_trip = original.vacuum_angular_wavenumber().vacuum_wavenumber();

        assert_eq!(
            round_trip.as_inverse_centimetres(),
            original.as_inverse_centimetres(),
        );
    }
}
