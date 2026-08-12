use std::fmt;

use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{LengthUnit, UnitLabel, VacuumAngularWavenumber, VacuumWavenumber};

/// Vacuum wavelength with an explicit physical length unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wavelength<F> {
    value: F,
    unit: LengthUnit,
}

impl<F> Wavelength<F> {
    /// Construct a vacuum wavelength.
    pub const fn new(value: F, unit: LengthUnit) -> Self {
        Self { value, unit }
    }

    pub fn metres(value: F) -> Self {
        Self::new(value, LengthUnit::Metre)
    }

    pub fn centimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Centimetre)
    }

    pub fn millimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Millimetre)
    }

    pub fn micrometres(value: F) -> Self {
        Self::new(value, LengthUnit::Micrometre)
    }

    pub fn nanometres(value: F) -> Self {
        Self::new(value, LengthUnit::Nanometre)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> LengthUnit {
        self.unit
    }

    pub fn into_parts(self) -> (F, LengthUnit) {
        (self.value, self.unit)
    }

    /// Return the wavelength expressed in centimetres.
    pub fn as_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_centimetres::<F>() * self.value
    }

    /// Convert to spectroscopic vacuum wavenumber.
    ///
    /// The relation is
    ///
    /// ```text
    /// ν̃ = 1 / λ₀.
    /// ```
    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FromPrimitive,
    {
        VacuumWavenumber::inverse_centimetres(F::one() / self.as_centimetres())
    }

    /// Convert to vacuum angular wavenumber.
    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FromPrimitive + FloatConst,
    {
        self.vacuum_wavenumber().vacuum_angular_wavenumber()
    }
}

impl<F> fmt::Display for Wavelength<F>
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
        let wavelength = Wavelength::nanometres(1550.0);

        assert_eq!(wavelength.value(), 1550.0);
        assert_eq!(wavelength.unit(), LengthUnit::Nanometre);
    }

    #[test]
    fn equivalent_units_have_equal_centimetre_values() {
        let wavelengths = [
            Wavelength::centimetres(1.5e-4),
            Wavelength::millimetres(1.5e-3),
            Wavelength::micrometres(1.5),
            Wavelength::nanometres(1500.0),
        ];

        for wavelength in wavelengths {
            approx::assert_relative_eq!(wavelength.as_centimetres(), 1.5e-4);
        }
    }

    #[test]
    fn into_parts_preserves_value_and_unit() {
        let wavelength = Wavelength::micrometres(2.0);

        assert_eq!(wavelength.into_parts(), (2.0, LengthUnit::Micrometre),);
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(Wavelength::nanometres(1550.0).to_string(), "1550 nm",);
    }

    #[test]
    fn wavelength_converts_to_vacuum_wavenumber() {
        let wavelength = Wavelength::micrometres(10.0);

        assert_eq!(
            wavelength.vacuum_wavenumber().as_inverse_centimetres(),
            1000.0,
        );
    }

    #[test]
    fn wavelength_converts_to_vacuum_angular_wavenumber() {
        let wavelength = Wavelength::micrometres(10.0);

        assert_eq!(
            wavelength
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
            std::f64::consts::TAU * 1000.0,
        );
    }

    #[test]
    fn wavelength_round_trips_through_vacuum_wavenumber() {
        let original = Wavelength::nanometres(1550.0);

        let round_trip = original.vacuum_wavenumber().wavelength();

        assert_eq!(round_trip.as_centimetres(), original.as_centimetres(),);
    }

    #[test]
    fn equivalent_spectral_quantities_produce_same_angular_wavenumber() {
        let wavelength = Wavelength::centimetres(1.0);

        let wavenumber = VacuumWavenumber::inverse_centimetres(1.0);

        assert_eq!(
            wavelength
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
            wavenumber
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
        );
    }
}
