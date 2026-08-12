use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{
    AngularFrequency, Frequency, VacuumAngularWavenumber, VacuumAngularWavenumberSquared,
    VacuumWavenumber, Wavelength,
};

/// A physical quantity that specifies a vacuum spectral coordinate.
///
/// Implementations convert their stored representation to vacuum angular
/// wavenumber,
///
/// ```text
/// k₀ = 2π / λ₀ = ω / c.
/// ```
///
/// The conversion preserves the physical meaning of the quantity while
/// expressing it in the spectral representation used throughout Lamina.
pub trait SpectralQuantity<F> {
    /// Convert this quantity to vacuum angular wavenumber.
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>;
}

impl<F> SpectralQuantity<F> for Wavelength<F>
where
    F: Float + FloatConst + FromPrimitive,
{
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F> {
        Wavelength::vacuum_angular_wavenumber(self)
    }
}

impl<F> SpectralQuantity<F> for VacuumWavenumber<F>
where
    F: Float + FloatConst + FromPrimitive,
{
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F> {
        VacuumWavenumber::vacuum_angular_wavenumber(self)
    }
}

impl<F> SpectralQuantity<F> for Frequency<F>
where
    F: Float + FloatConst + FromPrimitive,
{
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F> {
        Frequency::vacuum_angular_wavenumber(self)
    }
}

impl<F> SpectralQuantity<F> for AngularFrequency<F>
where
    F: Float + FloatConst + FromPrimitive,
{
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F> {
        AngularFrequency::vacuum_angular_wavenumber(self)
    }
}

impl<F> SpectralQuantity<F> for VacuumAngularWavenumber<F>
where
    F: Copy,
{
    fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F> {
        *self
    }
}

pub trait SpectralQuantityExt<F>: SpectralQuantity<F>
where
    F: Float + FloatConst + FromPrimitive,
{
    fn squared(&self) -> VacuumAngularWavenumberSquared<F> {
        let k0 = self.vacuum_angular_wavenumber().as_inverse_centimetres();

        VacuumAngularWavenumberSquared::inverse_centimetres_squared(k0 * k0)
    }
}

impl<F, T> SpectralQuantityExt<F> for T
where
    F: Float + FloatConst + FromPrimitive,
    T: SpectralQuantity<F>,
{
}

#[cfg(test)]
mod tests {
    use super::*;

    fn canonical<S>(spectral: &S) -> f64
    where
        S: SpectralQuantity<f64>,
    {
        spectral
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres()
    }

    #[test]
    fn equivalent_spectral_quantities_share_one_contract() {
        let expected = std::f64::consts::TAU;

        let wavelength = Wavelength::centimetres(1.0);
        let wavenumber = VacuumWavenumber::inverse_centimetres(1.0);
        let angular_wavenumber = VacuumAngularWavenumber::inverse_centimetres(expected);

        assert_eq!(canonical(&wavelength), expected);
        assert_eq!(canonical(&wavenumber), expected);
        assert_eq!(canonical(&angular_wavenumber), expected);
    }
}
