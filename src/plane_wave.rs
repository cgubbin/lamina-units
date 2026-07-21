use num_traits::{Float, FromPrimitive};

use crate::{
    AxisMetadata, RefractiveIndexJet,
    in_plane::{BivariateTransformJet, InPlaneCoordinate, InPlaneTransformError},
    spectral::{SpectralCoordinate, SpectralTransformError},
    transform::TransformJet,
};

/// A caller-facing pair of spectral and in-plane coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlaneWaveCoordinates {
    pub spectral: SpectralCoordinate,
    pub in_plane: InPlaneCoordinate,
}

/// Canonical plane-wave coordinates and their derivatives with respect to the
/// caller's spectral and in-plane coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneWaveCoordinateTransform<T> {
    /// Canonical vacuum angular wavenumber `k₀`.
    pub vacuum_angular_wavenumber: BivariateTransformJet<T>,

    /// Canonical parallel angular wavenumber `k∥`.
    pub parallel_angular_wavenumber: BivariateTransformJet<T>,
}

// #[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
// pub enum CoordinateTransformError {
//     #[error("failed to transform spectral coordinate: {0}")]
//     Spectral(#[from] SpectralTransformError),

//     #[error("failed to transform in-plane coordinate: {0}")]
//     InPlane(#[from] InPlaneTransformError),
// }

/// Errors produced while converting plane-wave coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CoordinateTransformError {
    #[error("coordinate value must be finite")]
    NonFinite,

    #[error("coordinate transformation is singular")]
    Singular,

    #[error("incident-angle coordinates require an incident refractive index")]
    MissingIncidentIndex,

    #[error("incident refractive index must be finite")]
    InvalidIncidentIndex,

    #[error("no real incident angle exists for these coordinates")]
    AngleOutsideDomain,
}

impl From<SpectralTransformError> for CoordinateTransformError {
    fn from(value: SpectralTransformError) -> Self {
        match value {
            SpectralTransformError::NonFinite => Self::NonFinite,
            SpectralTransformError::Singular => Self::Singular,
        }
    }
}

impl From<InPlaneTransformError> for CoordinateTransformError {
    fn from(value: InPlaneTransformError) -> Self {
        match value {
            InPlaneTransformError::NonFinite => Self::NonFinite,
            InPlaneTransformError::Singular => Self::Singular,
            InPlaneTransformError::MissingIncidentIndex => Self::MissingIncidentIndex,
            InPlaneTransformError::InvalidIncidentIndex => Self::InvalidIncidentIndex,
            InPlaneTransformError::AngleOutsideDomain => Self::AngleOutsideDomain,
        }
    }
}

impl PlaneWaveCoordinates {
    pub const fn new(spectral: SpectralCoordinate, in_plane: InPlaneCoordinate) -> Self {
        Self { spectral, in_plane }
    }

    pub const fn spectral(&self) -> SpectralCoordinate {
        self.spectral
    }

    pub const fn in_plane(&self) -> InPlaneCoordinate {
        self.in_plane
    }

    // pub fn requires_incident_index(&self) -> bool;

    pub fn transform_spectral<T>(
        &self,
        value: T,
    ) -> Result<TransformJet<T>, CoordinateTransformError>
    where
        T: Float + FromPrimitive,
    {
        Ok(self.spectral.to_canonical(value)?)
    }

    pub const fn requires_incident_index(&self) -> bool {
        matches!(self.in_plane, InPlaneCoordinate::IncidentAngle(_))
    }

    /// Converts caller-facing plane-wave coordinates into canonical backend
    /// coordinates.
    ///
    /// This method first transforms the spectral coordinate into canonical
    /// vacuum angular wavenumber and then resolves the in-plane coordinate
    /// using that spectral transformation.
    ///
    /// For incident-angle parameterisations, `context` must contain the
    /// incident refractive index and its first two derivatives with respect to
    /// canonical `k₀`.
    pub fn to_canonical<T>(
        self,
        spectral_value: T,
        in_plane_value: T,
    ) -> Result<PlaneWaveCoordinateTransform<T>, CoordinateTransformError>
    where
        T: Float + FromPrimitive,
    {
        let spectral = self.spectral.to_canonical(spectral_value)?;

        self.to_canonical_with_spectral(spectral, in_plane_value)
    }

    /// Converts caller-facing plane-wave coordinates into canonical backend
    /// coordinates.
    ///
    /// This method first transforms the spectral coordinate into canonical
    /// vacuum angular wavenumber and then resolves the in-plane coordinate
    /// using that spectral transformation.
    ///
    /// For incident-angle parameterisations, `context` must contain the
    /// incident refractive index and its first two derivatives with respect to
    /// canonical `k₀`.
    pub fn to_canonical_with_incident_index<T>(
        self,
        spectral: TransformJet<T>,
        in_plane_value: T,
        context: RefractiveIndexJet<T>,
    ) -> Result<PlaneWaveCoordinateTransform<T>, CoordinateTransformError>
    where
        T: Float + FromPrimitive,
    {
        self.resolve(spectral, in_plane_value, Some(context))
    }

    /// Converts plane-wave coordinates when the spectral transformation has
    /// already been evaluated.
    ///
    /// This is useful when the caller must evaluate a dispersive incident
    /// material at canonical `k₀` before resolving an incident angle.
    pub fn to_canonical_with_spectral<T>(
        self,
        spectral: TransformJet<T>,
        in_plane_value: T,
    ) -> Result<PlaneWaveCoordinateTransform<T>, CoordinateTransformError>
    where
        T: Float + FromPrimitive,
    {
        self.resolve(spectral, in_plane_value, None)
    }

    fn resolve<T>(
        self,
        spectral: TransformJet<T>,
        in_plane_value: T,
        context: Option<RefractiveIndexJet<T>>,
    ) -> Result<PlaneWaveCoordinateTransform<T>, CoordinateTransformError>
    where
        T: Float + FromPrimitive,
    {
        let zero = T::zero();

        let vacuum_angular_wavenumber = BivariateTransformJet::new(
            spectral.value,
            spectral.first,
            zero,
            spectral.second,
            zero,
            zero,
        );

        let parallel_angular_wavenumber =
            self.in_plane
                .to_canonical(in_plane_value, spectral, context)?;

        Ok(PlaneWaveCoordinateTransform {
            vacuum_angular_wavenumber,
            parallel_angular_wavenumber,
        })
    }

    pub fn axis_metadata(&self, axis: PlaneWaveAxis) -> AxisMetadata {
        match axis {
            PlaneWaveAxis::Spectral => self.spectral.metadata(),
            PlaneWaveAxis::InPlane => self.in_plane.metadata(),
        }
    }

    pub fn axes(&self) -> [AxisMetadata; 2] {
        [self.spectral.metadata(), self.in_plane.metadata()]
    }
}

/// Caller-facing axes of a plane-wave parameterisation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PlaneWaveAxis {
    Spectral,
    InPlane,
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::{
        angle::AngleUnit, frequency::FrequencyUnit, in_plane::InPlaneCoordinate,
        inverse_length::InverseLengthUnit, length::LengthUnit, metadata::AxisQuantity,
        spectral::SpectralCoordinate,
    };

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-13_f64.max(1.0e-12 * expected.abs());

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}; \
             absolute error was {:e}",
            (actual - expected).abs(),
        );
    }

    #[test]
    fn direct_parallel_angular_wavenumber_is_spectrally_independent() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
            InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
        );

        let transform = coordinates.to_canonical(10.0_f64, 500.0).unwrap();

        let parallel = transform.parallel_angular_wavenumber;

        assert_close(parallel.value, 500.0);
        assert_eq!(parallel.spectral_first, 0.0);
        assert_close(parallel.in_plane_first, 1.0);
        assert_eq!(parallel.spectral_second, 0.0);
        assert_eq!(parallel.mixed_second, 0.0);
        assert_eq!(parallel.in_plane_second, 0.0);
    }

    #[test]
    fn effective_index_couples_to_spectral_coordinate() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
            InPlaneCoordinate::EffectiveIndex,
        );

        let transform = coordinates.to_canonical(10.0_f64, 1.5).unwrap();

        let spectral = transform.vacuum_angular_wavenumber;
        let parallel = transform.parallel_angular_wavenumber;

        assert_close(parallel.value, 1.5 * spectral.value);
        assert_close(parallel.spectral_first, 1.5 * spectral.spectral_first);
        assert_close(parallel.in_plane_first, spectral.value);
        assert_close(parallel.spectral_second, 1.5 * spectral.spectral_second);
        assert_close(parallel.mixed_second, spectral.spectral_first);
        assert_eq!(parallel.in_plane_second, 0.0);
    }

    #[test]
    fn normal_incidence_has_zero_parallel_wavenumber() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
            InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
        );

        let index = TransformJet::new(1.5_f64, 0.0, 0.0);

        let spectral = coordinates.transform_spectral(30.0).unwrap();
        let transform = coordinates
            .to_canonical_with_incident_index(spectral, 0.0, index)
            .unwrap();

        let parallel = transform.parallel_angular_wavenumber;
        let spectral = transform.vacuum_angular_wavenumber;

        assert_close(parallel.value, 0.0);
        assert_close(parallel.spectral_first, 0.0);

        let radians_per_degree = std::f64::consts::PI / 180.0;

        assert_close(
            parallel.in_plane_first,
            1.5 * spectral.value * radians_per_degree,
        );
    }

    #[test]
    fn angle_in_degrees_and_radians_are_consistent() {
        let spectral = SpectralCoordinate::Frequency(FrequencyUnit::Terahertz);

        let degree_coordinates = PlaneWaveCoordinates::new(
            spectral,
            InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
        );

        let radian_coordinates = PlaneWaveCoordinates::new(
            spectral,
            InPlaneCoordinate::IncidentAngle(AngleUnit::Radian),
        );

        let context = TransformJet::new(1.5_f64, 0.0, 0.0);

        let spectral = degree_coordinates.transform_spectral(100.0).unwrap();
        let degree_transform = degree_coordinates
            .to_canonical_with_incident_index(spectral, 45.0, context)
            .unwrap();

        let spectral = radian_coordinates.transform_spectral(100.0).unwrap();
        let radian_transform = radian_coordinates
            .to_canonical_with_incident_index(spectral, std::f64::consts::FRAC_PI_4, context)
            .unwrap();

        assert_close(
            degree_transform.parallel_angular_wavenumber.value,
            radian_transform.parallel_angular_wavenumber.value,
        );
    }

    #[test]
    fn fixed_angle_spectral_derivative_includes_k0_dependence() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
            InPlaneCoordinate::IncidentAngle(AngleUnit::Radian),
        );

        let index = TransformJet::new(1.5_f64, 0.0, 0.0);

        let spectral = coordinates.transform_spectral(100.0).unwrap();
        let transform = coordinates
            .to_canonical_with_incident_index(spectral, std::f64::consts::FRAC_PI_6, index)
            .unwrap();

        let spectral = transform.vacuum_angular_wavenumber;
        let parallel = transform.parallel_angular_wavenumber;

        assert_close(parallel.spectral_first, 1.5 * spectral.spectral_first * 0.5);
    }

    #[test]
    fn incident_index_dispersion_contributes_to_spectral_derivative() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
            InPlaneCoordinate::IncidentAngle(AngleUnit::Radian),
        );

        let index = TransformJet::new(1.5_f64, 0.01, 0.002);

        let spectral = coordinates.transform_spectral(100.0).unwrap();

        let transform = coordinates
            .to_canonical_with_incident_index(spectral, std::f64::consts::FRAC_PI_6, index)
            .unwrap();

        let parallel = transform.parallel_angular_wavenumber;

        let expected_first = (1.5 + 100.0 * 0.01) * 0.5;

        let expected_second = (2.0 * 0.01 + 100.0 * 0.002) * 0.5;

        assert_close(parallel.spectral_first, expected_first);
        assert_close(parallel.spectral_second, expected_second);
    }

    #[test]
    fn incident_angle_requires_incident_index() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
            InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
        );

        let result = coordinates.to_canonical(30.0_f64, 45.0);

        assert!(matches!(
            result,
            Err(CoordinateTransformError::MissingIncidentIndex),
        ));
    }

    #[test]
    fn parallel_wavenumber_is_multiplied_by_two_pi() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
            InPlaneCoordinate::ParallelWavenumber(InverseLengthUnit::PerCentimetre),
        );

        let transform = coordinates.to_canonical(100.0_f64, 1.0).unwrap();

        let parallel = transform.parallel_angular_wavenumber;

        assert_close(parallel.value, std::f64::consts::TAU);
        assert_close(parallel.in_plane_first, std::f64::consts::TAU);
    }

    #[test]
    fn plane_wave_axes_report_caller_parameterisation() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
            InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
        );

        let [spectral, in_plane] = coordinates.axes();

        assert_eq!(spectral.quantity, AxisQuantity::VacuumWavelength,);
        assert_eq!(spectral.quantity_name, "Vacuum wavelength",);
        assert_eq!(spectral.quantity_symbol, "λ₀");
        assert_eq!(spectral.unit_symbol, Some("µm"));
        assert_eq!(spectral.label, "λ₀ (µm)");

        assert_eq!(in_plane.quantity, AxisQuantity::IncidentAngle,);
        assert_eq!(in_plane.quantity_name, "Incident angle",);
        assert_eq!(in_plane.quantity_symbol, "θ");
        assert_eq!(in_plane.unit_symbol, Some("°"));
        assert_eq!(in_plane.label, "θ (°)");
    }

    #[test]
    fn effective_index_axis_has_no_unit() {
        let coordinates = PlaneWaveCoordinates::new(
            SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
            InPlaneCoordinate::EffectiveIndex,
        );

        let metadata = coordinates.axis_metadata(PlaneWaveAxis::InPlane);

        assert_eq!(metadata.quantity, AxisQuantity::EffectiveIndex,);
        assert_eq!(metadata.unit_symbol, None);
        assert_eq!(metadata.label, "n_eff");
    }

    #[test]
    fn axis_metadata_preserves_angular_and_cyclic_distinction() {
        let angular = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
            InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
        );

        let cyclic = PlaneWaveCoordinates::new(
            SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre),
            InPlaneCoordinate::ParallelWavenumber(InverseLengthUnit::PerCentimetre),
        );

        assert_eq!(angular.axes()[0].label, "k₀ (cm⁻¹)",);
        assert_eq!(angular.axes()[1].label, "k∥ (cm⁻¹)",);

        assert_eq!(cyclic.axes()[0].label, "ν̃ (cm⁻¹)",);
        assert_eq!(cyclic.axes()[1].label, "ν̃∥ (cm⁻¹)",);
    }
}
