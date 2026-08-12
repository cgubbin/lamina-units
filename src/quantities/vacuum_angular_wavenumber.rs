use std::fmt;

use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{InverseLengthUnit, UnitLabel, VacuumWavenumber, Wavelength};

/// Vacuum angular wavenumber.
///
/// The vacuum angular wavenumber is
///
/// ```text
/// k₀ = 2π / λ₀ = ω / c.
/// ```
///
/// It has dimensions of inverse length. Lamina uses inverse centimetres as
/// the reference representation when a numerical `k₀` is required.
///
/// This quantity is distinct from [`VacuumWavenumber`](crate::VacuumWavenumber),
/// which represents the spectroscopic wavenumber `1 / λ₀`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VacuumAngularWavenumber<F> {
    value: F,
    unit: InverseLengthUnit,
}

impl<F> VacuumAngularWavenumber<F> {
    /// Construct a vacuum angular wavenumber.
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

    /// Return the numerical value in the stored unit.
    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    /// Return the stored inverse-length unit.
    pub const fn unit(&self) -> InverseLengthUnit {
        self.unit
    }

    /// Consume this quantity and return its value and unit.
    pub fn into_parts(self) -> (F, InverseLengthUnit) {
        (self.value, self.unit)
    }

    /// Return the vacuum angular wavenumber in inverse centimetres.
    pub fn as_inverse_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_inverse_centimetres::<F>() * self.value
    }

    /// Convert to spectroscopic vacuum wavenumber.
    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        VacuumWavenumber::inverse_centimetres(
            self.as_inverse_centimetres() / crate::private::two_pi(),
        )
    }

    /// Convert to vacuum wavelength.
    pub fn wavelength(&self) -> Wavelength<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        self.vacuum_wavenumber().wavelength()
    }
}

impl<F> fmt::Display for VacuumAngularWavenumber<F>
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
        let k0 = VacuumAngularWavenumber::inverse_micrometres(2.5);

        assert_eq!(k0.value(), 2.5);
        assert_eq!(k0.unit(), InverseLengthUnit::PerMicrometre);
    }

    #[test]
    fn equivalent_units_have_equal_inverse_centimetre_values() {
        let values = [
            VacuumAngularWavenumber::inverse_metres(1.0e5),
            VacuumAngularWavenumber::inverse_centimetres(1.0e3),
            VacuumAngularWavenumber::inverse_millimetres(100.0),
            VacuumAngularWavenumber::inverse_micrometres(0.1),
            VacuumAngularWavenumber::inverse_nanometres(1.0e-4),
        ];

        for value in values {
            assert_eq!(value.as_inverse_centimetres(), 1.0e3);
        }
    }

    #[test]
    fn into_parts_preserves_value_and_unit() {
        let k0 = VacuumAngularWavenumber::inverse_centimetres(2500.0);

        assert_eq!(k0.into_parts(), (2500.0, InverseLengthUnit::PerCentimetre),);
    }

    #[test]
    fn display_includes_value_and_unit() {
        assert_eq!(
            VacuumAngularWavenumber::inverse_centimetres(2500.0).to_string(),
            "2500 cm⁻¹",
        );
    }
}
