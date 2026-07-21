use core::option::Option::None;

use crate::{
    AxisMetadata,
    angle::AngleUnit,
    inverse_length::InverseLengthUnit,
    label::{CoordinateLabel, UnitLabel},
    metadata::AxisQuantity,
    private::conversion_constant,
    transform::{TransformJet, UnitScale},
};

use num_traits::{Float, FromPrimitive};

/// The caller-facing parameterisation of the in-plane wavevector.
///
/// These variants ultimately map to the canonical parallel angular
/// wavenumber `k∥`, expressed in inverse centimetres.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum InPlaneCoordinate {
    /// Parallel angular wavenumber `k∥`.
    ParallelAngularWavenumber(InverseLengthUnit),

    /// Parallel cyclic wavenumber `k∥ / 2π`.
    ParallelWavenumber(InverseLengthUnit),

    /// Effective index `n_eff = k∥ / k₀`.
    EffectiveIndex,

    /// Angle of incidence in the incident medium.
    IncidentAngle(AngleUnit),
}

/// Refractive index and its first two derivatives with respect to canonical
/// vacuum angular wavenumber `k₀`.
///
/// The fields of the underlying [`TransformJet`] have the meaning
///
/// ```text
/// value  = n(k₀)
/// first  = dn/dk₀
/// second = d²n/dk₀²
/// ```
pub type RefractiveIndexJet<T> = TransformJet<T>;

/// The value and derivatives of a scalar canonical coordinate with respect to
/// the caller's spectral and in-plane coordinates.
///
/// If the caller coordinates are `(s, p)`, the fields represent
///
/// ```text
/// value             = q(s, p)
/// spectral_first    = ∂q/∂s
/// in_plane_first    = ∂q/∂p
/// spectral_second   = ∂²q/∂s²
/// mixed_second      = ∂²q/(∂s ∂p)
/// in_plane_second   = ∂²q/∂p²
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BivariateTransformJet<T> {
    pub value: T,
    pub spectral_first: T,
    pub in_plane_first: T,
    pub spectral_second: T,
    pub mixed_second: T,
    pub in_plane_second: T,
}

impl<T> BivariateTransformJet<T> {
    pub const fn new(
        value: T,
        spectral_first: T,
        in_plane_first: T,
        spectral_second: T,
        mixed_second: T,
        in_plane_second: T,
    ) -> Self {
        Self {
            value,
            spectral_first,
            in_plane_first,
            spectral_second,
            mixed_second,
            in_plane_second,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum InPlaneTransformError {
    #[error("the in-plane coordinate transformation requires finite values")]
    NonFinite,

    #[error("incident-angle parameterisation requires an incident refractive index")]
    MissingIncidentIndex,

    #[error("incident refractive index must be finite")]
    InvalidIncidentIndex,

    #[error("the incident-angle transformation is singular")]
    Singular,

    #[error("no real principal incident angle exists")]
    AngleOutsideDomain,
}

impl InPlaneCoordinate {
    /// Transforms the caller's in-plane coordinate into canonical parallel
    /// angular wavenumber.
    ///
    /// `spectral` must contain the canonical vacuum angular wavenumber and its
    /// derivatives with respect to the caller's spectral coordinate:
    ///
    /// ```text
    /// value  = k₀(s)
    /// first  = dk₀/ds
    /// second = d²k₀/ds²
    /// ```
    ///
    /// The returned bivariate jet describes `k∥(s, p)`, where `p` is the
    /// caller's in-plane coordinate.
    pub fn to_canonical<T>(
        self,
        in_plane_value: T,
        spectral: TransformJet<T>,
        context: Option<RefractiveIndexJet<T>>,
    ) -> Result<BivariateTransformJet<T>, InPlaneTransformError>
    where
        T: Float + FromPrimitive,
    {
        if !in_plane_value.is_finite()
            || !spectral.value.is_finite()
            || !spectral.first.is_finite()
            || !spectral.second.is_finite()
        {
            return Err(InPlaneTransformError::NonFinite);
        }

        match self {
            Self::ParallelAngularWavenumber(unit) => {
                let scale = unit.scale();

                Ok(BivariateTransformJet::new(
                    in_plane_value * scale,
                    T::zero(),
                    scale,
                    T::zero(),
                    T::zero(),
                    T::zero(),
                ))
            }

            Self::ParallelWavenumber(unit) => {
                let tau: T = conversion_constant(std::f64::consts::TAU);
                let scale = tau * unit.scale();

                Ok(BivariateTransformJet::new(
                    in_plane_value * scale,
                    T::zero(),
                    scale,
                    T::zero(),
                    T::zero(),
                    T::zero(),
                ))
            }

            Self::EffectiveIndex => {
                let effective_index = in_plane_value;
                let k0 = spectral.value;
                let k0_first = spectral.first;
                let k0_second = spectral.second;

                Ok(BivariateTransformJet::new(
                    effective_index * k0,
                    effective_index * k0_first,
                    k0,
                    effective_index * k0_second,
                    k0_first,
                    T::zero(),
                ))
            }

            Self::IncidentAngle(unit) => {
                let Some(index) = context else {
                    return Err(InPlaneTransformError::MissingIncidentIndex);
                };

                if !index.value.is_finite() || !index.first.is_finite() || !index.second.is_finite()
                {
                    return Err(InPlaneTransformError::InvalidIncidentIndex);
                }

                let angle_scale = unit.scale();
                let angle = in_plane_value * angle_scale;
                let sin_angle = angle.sin();
                let cos_angle = angle.cos();

                let k0 = spectral.value;
                let k0_first = spectral.first;
                let k0_second = spectral.second;

                let index_value = index.value;
                let index_first = index.first;
                let index_second = index.second;

                // h(k₀) = k₀ n(k₀)
                let h = k0 * index_value;
                let h_first = index_value + k0 * index_first;
                let h_second = (T::one() + T::one()) * index_first + k0 * index_second;

                let value = h * sin_angle;

                let spectral_first = h_first * k0_first * sin_angle;

                let in_plane_first = h * cos_angle * angle_scale;

                let spectral_second =
                    (h_second * k0_first * k0_first + h_first * k0_second) * sin_angle;

                let mixed_second = h_first * k0_first * cos_angle * angle_scale;

                let in_plane_second = -h * sin_angle * angle_scale * angle_scale;

                Ok(BivariateTransformJet::new(
                    value,
                    spectral_first,
                    in_plane_first,
                    spectral_second,
                    mixed_second,
                    in_plane_second,
                ))
            }
        }
    }

    /// Converts canonical `k∥` back to the selected caller coordinate.
    ///
    /// The inverse mapping is directly available for:
    ///
    /// - parallel angular wavenumber;
    /// - parallel cyclic wavenumber;
    /// - effective index.
    ///
    /// Incident angle additionally requires the incident refractive index.
    ///
    /// The returned jet contains derivatives with respect to canonical
    /// parallel angular wavenumber while holding `k₀` and the supplied
    /// incident index fixed.
    pub fn from_canonical<T>(
        self,
        parallel_angular_wavenumber: T,
        vacuum_angular_wavenumber: T,
        context: Option<RefractiveIndexJet<T>>,
    ) -> Result<TransformJet<T>, InPlaneTransformError>
    where
        T: Float + FromPrimitive,
    {
        if !parallel_angular_wavenumber.is_finite() || !vacuum_angular_wavenumber.is_finite() {
            return Err(InPlaneTransformError::NonFinite);
        }

        match self {
            Self::ParallelAngularWavenumber(unit) => {
                let inverse_scale = T::one() / unit.scale();

                Ok(TransformJet::new(
                    parallel_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::ParallelWavenumber(unit) => {
                let tau: T = conversion_constant(std::f64::consts::TAU);
                let inverse_scale = T::one() / (tau * unit.scale());

                Ok(TransformJet::new(
                    parallel_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::EffectiveIndex => {
                let inverse_scale = T::one() / vacuum_angular_wavenumber;

                Ok(TransformJet::new(
                    parallel_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::IncidentAngle(unit) => {
                let Some(index) = context else {
                    return Err(InPlaneTransformError::MissingIncidentIndex);
                };

                if !index.value.is_finite() {
                    return Err(InPlaneTransformError::InvalidIncidentIndex);
                }

                let denominator = vacuum_angular_wavenumber * index.value;

                if denominator == T::zero() {
                    return Err(InPlaneTransformError::Singular);
                }

                let normalised = parallel_angular_wavenumber / denominator;

                if normalised < -T::one() || normalised > T::one() {
                    return Err(InPlaneTransformError::AngleOutsideDomain);
                }

                let angle_radians = normalised.asin();
                let root = (T::one() - normalised * normalised).sqrt();

                if root == T::zero() {
                    return Err(InPlaneTransformError::Singular);
                }

                let inverse_unit_scale = T::one() / unit.scale();

                let first = inverse_unit_scale / (denominator * root);

                let second = inverse_unit_scale * normalised
                    / (denominator * denominator * root * root * root);

                Ok(TransformJet::new(
                    angle_radians * inverse_unit_scale,
                    first,
                    second,
                ))
            }
        }
    }

    pub(crate) fn metadata(&self) -> AxisMetadata {
        AxisMetadata {
            quantity: match self {
                Self::ParallelAngularWavenumber(_) => AxisQuantity::ParallelAngularWavenumber,
                Self::ParallelWavenumber(_) => AxisQuantity::ParallelWavenumber,
                Self::EffectiveIndex => AxisQuantity::EffectiveIndex,
                Self::IncidentAngle(_) => AxisQuantity::IncidentAngle,
            },
            quantity_name: self.quantity_name(),
            quantity_symbol: self.quantity_symbol(),
            unit_symbol: {
                let symbol = self.unit_symbol();
                if symbol.is_empty() {
                    None
                } else {
                    Some(symbol)
                }
            },
            label: self.axis_label(),
        }
    }
}

impl CoordinateLabel for InPlaneCoordinate {
    fn quantity_name(&self) -> &'static str {
        match self {
            Self::ParallelAngularWavenumber(_) => "Parallel angular wavenumber",
            Self::ParallelWavenumber(_) => "Parallel wavenumber",
            Self::EffectiveIndex => "Effective index",
            Self::IncidentAngle(_) => "Incident angle",
        }
    }

    fn quantity_symbol(&self) -> &'static str {
        match self {
            Self::ParallelAngularWavenumber(_) => "k∥",
            Self::ParallelWavenumber(_) => "ν̃∥",
            Self::EffectiveIndex => "n_eff",
            Self::IncidentAngle(_) => "θ",
        }
    }

    fn unit_symbol(&self) -> &'static str {
        match self {
            Self::ParallelAngularWavenumber(unit) | Self::ParallelWavenumber(unit) => unit.symbol(),
            Self::EffectiveIndex => "",
            Self::IncidentAngle(unit) => unit.symbol(),
        }
    }

    fn axis_label(&self) -> String {
        match self {
            Self::EffectiveIndex => self.quantity_symbol().to_owned(),
            _ => format!("{} ({})", self.quantity_symbol(), self.unit_symbol(),),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::{angle::AngleUnit, inverse_length::InverseLengthUnit, transform::TransformJet};

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-13_f64.max(expected.abs() * 1.0e-12);

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}; \
             error was {:e}",
            (actual - expected).abs(),
        );
    }

    fn identity_spectral(k0: f64) -> TransformJet<f64> {
        TransformJet::new(k0, 1.0, 0.0)
    }

    #[test]
    fn parallel_angular_wavenumber_in_inverse_centimetres_is_identity() {
        let coordinate =
            InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre);

        let jet = coordinate
            .to_canonical(12.5_f64, identity_spectral(100.0), None)
            .unwrap();

        assert_close(jet.value, 12.5);
        assert_eq!(jet.spectral_first, 0.0);
        assert_close(jet.in_plane_first, 1.0);
        assert_eq!(jet.spectral_second, 0.0);
        assert_eq!(jet.mixed_second, 0.0);
        assert_eq!(jet.in_plane_second, 0.0);
    }

    #[test]
    fn parallel_angular_wavenumber_applies_unit_scale() {
        let coordinate = InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerMetre);

        let jet = coordinate
            .to_canonical(500.0_f64, identity_spectral(100.0), None)
            .unwrap();

        assert_close(jet.value, 5.0);
        assert_close(jet.in_plane_first, 1.0e-2);
    }

    #[test]
    fn parallel_wavenumber_includes_two_pi() {
        let coordinate = InPlaneCoordinate::ParallelWavenumber(InverseLengthUnit::PerCentimetre);

        let jet = coordinate
            .to_canonical(1.0_f64, identity_spectral(100.0), None)
            .unwrap();

        assert_close(jet.value, std::f64::consts::TAU);

        assert_close(jet.in_plane_first, std::f64::consts::TAU);
    }

    #[test]
    fn direct_parallel_coordinates_do_not_depend_on_spectral_coordinate() {
        let coordinate =
            InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre);

        let first = coordinate
            .to_canonical(12.0_f64, TransformJet::new(10.0, 2.0, 3.0), None)
            .unwrap();

        let second = coordinate
            .to_canonical(12.0_f64, TransformJet::new(100.0, -7.0, 9.0), None)
            .unwrap();

        assert_eq!(first, second);
    }

    #[test]
    fn effective_index_value_is_neff_times_k0() {
        let spectral = TransformJet::new(100.0_f64, 2.0, 3.0);

        let jet = InPlaneCoordinate::EffectiveIndex
            .to_canonical(1.5, spectral, None)
            .unwrap();

        assert_close(jet.value, 150.0);
    }

    #[test]
    fn effective_index_derivatives_are_correct() {
        let spectral = TransformJet::new(100.0_f64, 2.0, 3.0);

        let jet = InPlaneCoordinate::EffectiveIndex
            .to_canonical(1.5, spectral, None)
            .unwrap();

        assert_close(jet.spectral_first, 3.0);
        assert_close(jet.in_plane_first, 100.0);

        assert_close(jet.spectral_second, 4.5);
        assert_close(jet.mixed_second, 2.0);
        assert_eq!(jet.in_plane_second, 0.0);
    }

    #[test]
    fn normal_incidence_has_zero_parallel_wavenumber() {
        let index = TransformJet::new(1.5_f64, 0.0, 0.0);

        let jet = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(0.0, identity_spectral(100.0), Some(index))
            .unwrap();

        assert_close(jet.value, 0.0);
        assert_close(jet.spectral_first, 0.0);
        assert_close(jet.spectral_second, 0.0);
        assert_close(jet.in_plane_second, 0.0);

        assert_close(jet.in_plane_first, 150.0);
        assert_close(jet.mixed_second, 1.5);
    }

    #[test]
    fn angle_value_is_n_k0_sin_theta() {
        let index = TransformJet::new(1.5_f64, 0.0, 0.0);

        let jet = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(
                std::f64::consts::FRAC_PI_6,
                identity_spectral(100.0),
                Some(index),
            )
            .unwrap();

        assert_close(jet.value, 75.0);
    }

    #[test]
    fn degree_and_radian_angles_produce_same_value() {
        let spectral = identity_spectral(100.0);

        let context = Some(TransformJet::new(1.5_f64, 0.0, 0.0));

        let degrees = InPlaneCoordinate::IncidentAngle(AngleUnit::Degree)
            .to_canonical(30.0, spectral, context)
            .unwrap();

        let radians = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(std::f64::consts::FRAC_PI_6, spectral, context)
            .unwrap();

        assert_close(degrees.value, radians.value);
        assert_close(degrees.spectral_first, radians.spectral_first);
        assert_close(degrees.spectral_second, radians.spectral_second);
    }

    #[test]
    fn angular_derivatives_include_angle_unit_scale() {
        let spectral = identity_spectral(100.0);

        let context = Some(TransformJet::new(1.5_f64, 0.0, 0.0));

        let degrees = InPlaneCoordinate::IncidentAngle(AngleUnit::Degree)
            .to_canonical(30.0, spectral, context)
            .unwrap();

        let radians = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(std::f64::consts::FRAC_PI_6, spectral, context)
            .unwrap();

        let radians_per_degree = std::f64::consts::PI / 180.0;

        assert_close(
            degrees.in_plane_first,
            radians.in_plane_first * radians_per_degree,
        );

        assert_close(
            degrees.mixed_second,
            radians.mixed_second * radians_per_degree,
        );

        assert_close(
            degrees.in_plane_second,
            radians.in_plane_second * radians_per_degree.powi(2),
        );
    }

    #[test]
    fn incident_index_dispersion_contributes_to_spectral_derivatives() {
        let spectral = TransformJet::new(100.0_f64, 2.0, 3.0);

        let index = TransformJet::new(1.5_f64, 0.01, 0.002);

        let angle = std::f64::consts::FRAC_PI_6;

        let jet = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(angle, spectral, Some(index))
            .unwrap();

        let sin_angle = angle.sin();

        let h_first = 1.5 + 100.0 * 0.01;

        let h_second = 2.0 * 0.01 + 100.0 * 0.002;

        let expected_first = h_first * 2.0 * sin_angle;

        let expected_second = (h_second * 2.0_f64.powi(2) + h_first * 3.0) * sin_angle;

        assert_close(jet.spectral_first, expected_first);

        assert_close(jet.spectral_second, expected_second);
    }

    #[test]
    fn mixed_angle_derivative_is_correct() {
        let spectral = TransformJet::new(100.0_f64, 2.0, 0.0);

        let index = TransformJet::new(1.5_f64, 0.01, 0.0);

        let angle = std::f64::consts::FRAC_PI_3;

        let jet = InPlaneCoordinate::IncidentAngle(AngleUnit::Radian)
            .to_canonical(angle, spectral, Some(index))
            .unwrap();

        let h_first = 1.5 + 100.0 * 0.01;

        let expected = h_first * 2.0 * angle.cos();

        assert_close(jet.mixed_second, expected);
    }

    #[test]
    fn incident_angle_requires_incident_index() {
        let result = InPlaneCoordinate::IncidentAngle(AngleUnit::Degree).to_canonical(
            30.0_f64,
            identity_spectral(100.0),
            None,
        );

        assert_eq!(result, Err(InPlaneTransformError::MissingIncidentIndex),);
    }

    #[test]
    fn non_finite_in_plane_value_is_rejected() {
        let result = InPlaneCoordinate::EffectiveIndex.to_canonical(
            f64::NAN,
            identity_spectral(100.0),
            None,
        );

        assert_eq!(result, Err(InPlaneTransformError::NonFinite),);
    }

    #[test]
    fn non_finite_spectral_jet_is_rejected() {
        let result = InPlaneCoordinate::EffectiveIndex.to_canonical(
            1.5_f64,
            TransformJet::new(100.0, f64::INFINITY, 0.0),
            None,
        );

        assert_eq!(result, Err(InPlaneTransformError::NonFinite),);
    }

    #[test]
    fn non_finite_incident_index_jet_is_rejected() {
        let result = InPlaneCoordinate::IncidentAngle(AngleUnit::Degree).to_canonical(
            30.0_f64,
            identity_spectral(100.0),
            Some(TransformJet::new(1.5, f64::NAN, 0.0)),
        );

        assert_eq!(result, Err(InPlaneTransformError::InvalidIncidentIndex),);
    }

    #[test]
    fn in_plane_transform_supports_f32() {
        let jet = InPlaneCoordinate::EffectiveIndex
            .to_canonical(1.5_f32, TransformJet::new(100.0_f32, 2.0, 3.0), None)
            .unwrap();

        assert_eq!(jet.value, 150.0);
        assert_eq!(jet.spectral_first, 3.0);
        assert_eq!(jet.in_plane_first, 100.0);
    }
}
