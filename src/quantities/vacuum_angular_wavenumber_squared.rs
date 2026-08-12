/// Squared vacuum angular wavenumber `k₀²`.
///
/// Lamina represents this quantity in inverse centimetres squared (`cm⁻²`).
///
/// This is useful for parameters that enter constitutive models naturally
/// as a squared spectral scale, such as oscillator strengths.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VacuumAngularWavenumberSquared<F> {
    value_cm2: F,
}

impl<F> VacuumAngularWavenumberSquared<F> {
    /// Construct directly from `k₀²` expressed in `cm⁻²`.
    pub const fn inverse_centimetres_squared(value: F) -> Self {
        Self { value_cm2: value }
    }

    /// Return `k₀²` expressed in `cm⁻²`.
    pub fn as_inverse_centimetres_squared(self) -> F {
        self.value_cm2
    }
}
