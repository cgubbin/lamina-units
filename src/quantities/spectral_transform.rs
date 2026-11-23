use num_traits::{Float, FloatConst, FromPrimitive};

use crate::{AngularFrequencyUnit, FrequencyUnit, InverseLengthUnit, LengthUnit};

use crate::private::two_pi;

use crate::constants::SPEED_OF_LIGHT_CM_PER_SECOND;

/// Intrinsic transformation from a caller-facing spectral coordinate to
/// vacuum angular wavenumber.
///
/// This type contains only the physical scalar transformation. It carries no
/// derivative or jet information; higher-level crates may apply the same
/// transformation using their own scalar or automatic-differentiation
/// algebra.
///
/// The canonical target quantity is vacuum angular wavenumber
///
/// ```text
/// k₀ = 2π / λ₀ = ω / c,
/// ```
///
/// expressed in inverse centimetres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectralTransform<R> {
    /// Linear transformation
    ///
    /// ```text
    /// k₀ = scale * x.
    /// ```
    Linear { scale: R },

    /// Reciprocal transformation
    ///
    /// ```text
    /// k₀ = numerator / (input_scale * x).
    /// ```
    ///
    /// `input_scale` converts the supplied quantity into the reference unit
    /// required by the reciprocal relation.
    Reciprocal { input_scale: R, numerator: R },
}

/// Spectral coordinate representation.
///
/// This identifies the physical quantity and unit used by a caller to specify
/// a vacuum spectral coordinate.
///
/// Conversion to vacuum angular wavenumber is intrinsic and does not depend on
/// an optical stack, material model, incidence direction, or other
/// problem-specific context.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpectralCoordinate {
    /// Vacuum angular wavenumber
    ///
    /// ```text
    /// k₀ = 2π / λ₀.
    /// ```
    VacuumAngularWavenumber(InverseLengthUnit),

    /// Spectroscopic vacuum wavenumber
    ///
    /// ```text
    /// ν̃ = 1 / λ₀.
    /// ```
    VacuumWavenumber(InverseLengthUnit),

    /// Ordinary frequency `f`.
    Frequency(FrequencyUnit),

    /// Angular frequency
    ///
    /// ```text
    /// ω = 2πf.
    /// ```
    AngularFrequency(AngularFrequencyUnit),

    /// Vacuum wavelength `λ₀`.
    VacuumWavelength(LengthUnit),
}

impl SpectralCoordinate {
    /// Return the intrinsic transformation from this coordinate to vacuum
    /// angular wavenumber in inverse centimetres.
    pub fn transform<R>(self) -> SpectralTransform<R>
    where
        R: Float + FloatConst + FromPrimitive,
    {
        let two_pi = two_pi::<R>();

        match self {
            Self::VacuumAngularWavenumber(unit) => SpectralTransform::Linear {
                scale: unit.scale_to_inverse_centimetres(),
            },

            Self::VacuumWavenumber(unit) => SpectralTransform::Linear {
                scale: two_pi * unit.scale_to_inverse_centimetres(),
            },

            Self::Frequency(unit) => {
                let c = R::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND)
                    .expect("speed of light must be representable");

                SpectralTransform::Linear {
                    scale: two_pi * unit.scale_to_hertz() / c,
                }
            }

            Self::AngularFrequency(unit) => {
                let c = R::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND)
                    .expect("speed of light must be representable");

                SpectralTransform::Linear {
                    scale: unit.scale_to_radians_per_second::<R>() / c,
                }
            }

            Self::VacuumWavelength(unit) => SpectralTransform::Reciprocal {
                input_scale: unit.scale_to_centimetres(),
                numerator: two_pi,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    const TOLERANCE: f64 = 1.0e-12;

    fn apply(transform: SpectralTransform<f64>, value: f64) -> f64 {
        match transform {
            SpectralTransform::Linear { scale } => scale * value,

            SpectralTransform::Reciprocal {
                input_scale,
                numerator,
            } => numerator / (input_scale * value),
        }
    }

    #[test]
    fn vacuum_angular_wavenumber_is_only_unit_scaled() {
        let coordinate =
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerMicrometre);

        let actual = apply(coordinate.transform(), 2.0);

        let expected = 2.0e4;

        assert_relative_eq!(
            actual,
            expected,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn vacuum_wavenumber_acquires_two_pi_factor() {
        let coordinate = SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre);

        let actual = apply(coordinate.transform(), 1000.0);

        let expected = std::f64::consts::TAU * 1000.0;

        assert_relative_eq!(
            actual,
            expected,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn frequency_uses_two_pi_over_speed_of_light() {
        let coordinate = SpectralCoordinate::Frequency(FrequencyUnit::Hertz);

        let actual = apply(coordinate.transform(), SPEED_OF_LIGHT_CM_PER_SECOND);

        assert_relative_eq!(
            actual,
            std::f64::consts::TAU,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn angular_frequency_uses_inverse_speed_of_light() {
        let coordinate =
            SpectralCoordinate::AngularFrequency(AngularFrequencyUnit::RadiansPerSecond);

        let actual = apply(
            coordinate.transform(),
            std::f64::consts::TAU * SPEED_OF_LIGHT_CM_PER_SECOND,
        );

        assert_relative_eq!(
            actual,
            std::f64::consts::TAU,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn wavelength_uses_reciprocal_relation() {
        let coordinate = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        let actual = apply(coordinate.transform(), 10.0);

        let expected = std::f64::consts::TAU / 1.0e-3;

        assert_relative_eq!(
            actual,
            expected,
            epsilon = TOLERANCE,
            max_relative = TOLERANCE,
        );
    }

    #[test]
    fn equivalent_spectral_coordinates_give_same_transform_value() {
        let wavelength_um = 10.0;

        let vacuum_wavenumber_cm = 1000.0;

        let frequency_hz = SPEED_OF_LIGHT_CM_PER_SECOND * vacuum_wavenumber_cm;

        let angular_frequency = std::f64::consts::TAU * frequency_hz;

        let vacuum_angular_wavenumber = std::f64::consts::TAU * vacuum_wavenumber_cm;

        let values = [
            apply(
                SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre).transform(),
                wavelength_um,
            ),
            apply(
                SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre).transform(),
                vacuum_wavenumber_cm,
            ),
            apply(
                SpectralCoordinate::Frequency(FrequencyUnit::Hertz).transform(),
                frequency_hz,
            ),
            apply(
                SpectralCoordinate::AngularFrequency(AngularFrequencyUnit::RadiansPerSecond)
                    .transform(),
                angular_frequency,
            ),
            apply(
                SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre)
                    .transform(),
                vacuum_angular_wavenumber,
            ),
        ];

        let expected = std::f64::consts::TAU * vacuum_wavenumber_cm;

        for actual in values {
            assert_relative_eq!(
                actual,
                expected,
                epsilon = TOLERANCE,
                max_relative = TOLERANCE,
            );
        }
    }
}
