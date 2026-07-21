use num_traits::{Float, FromPrimitive};

use crate::{
    AxisMetadata,
    frequency::{AngularFrequencyUnit, FrequencyUnit, SPEED_OF_LIGHT_CM_PER_S},
    inverse_length::InverseLengthUnit,
    label::{CoordinateLabel, UnitLabel},
    length::LengthUnit,
    metadata::AxisQuantity,
    private::conversion_constant,
    transform::{TransformJet, UnitScale},
};

/// A spectral coordinate used to parameterise an optical calculation.
///
/// Every variant is transformed into the canonical coordinate expected by the
/// TMM backend: the vacuum angular wavenumber
///
/// ```text
/// k₀ = 2π / λ,
/// ```
///
/// expressed in inverse centimetres.
///
/// [`VacuumWavenumber`](Self::VacuumWavenumber) denotes the spectroscopic
/// vacuum wavenumber
///
/// ```text
/// ν̃ = 1 / λ,
/// ```
///
/// and is therefore distinct from
/// [`VacuumAngularWavenumber`](Self::VacuumAngularWavenumber).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SpectralCoordinate {
    /// Vacuum angular wavenumber `k₀ = 2π / λ`.
    VacuumAngularWavenumber(InverseLengthUnit),

    /// Spectroscopic vacuum wavenumber `ν̃ = 1 / λ`.
    VacuumWavenumber(InverseLengthUnit),

    /// Ordinary frequency `f`.
    Frequency(FrequencyUnit),

    /// Angular frequency `ω = 2πf`.
    AngularFrequency(AngularFrequencyUnit),

    /// Vacuum wavelength `λ`.
    VacuumWavelength(LengthUnit),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum SpectralTransformError {
    #[error("coordinate value must be finite")]
    NonFinite,

    #[error("coordinate value must not be singular")]
    Singular,
}

impl SpectralCoordinate {
    pub(crate) fn to_canonical<T>(self, value: T) -> Result<TransformJet<T>, SpectralTransformError>
    where
        T: Float + FromPrimitive,
    {
        if !value.is_finite() {
            return Err(SpectralTransformError::NonFinite);
        }

        let tau: T = conversion_constant(std::f64::consts::TAU);
        let c: T = conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

        match self {
            Self::VacuumAngularWavenumber(unit) => {
                let scale = unit.scale();

                Ok(TransformJet::new(value * scale, scale, T::zero()))
            }

            Self::VacuumWavenumber(unit) => {
                let scale = tau * unit.scale();

                Ok(TransformJet::new(value * scale, scale, T::zero()))
            }

            Self::Frequency(unit) => {
                let scale = tau * unit.scale() / c;

                Ok(TransformJet::new(value * scale, scale, T::zero()))
            }

            Self::AngularFrequency(unit) => {
                let scale: T = <AngularFrequencyUnit as UnitScale<T>>::scale(unit) / c;

                Ok(TransformJet::new(value * scale, scale, T::zero()))
            }

            Self::VacuumWavelength(unit) => {
                let length_scale = unit.scale();
                let wavelength_cm = value * length_scale;

                if !wavelength_cm.is_finite() {
                    return Err(SpectralTransformError::NonFinite);
                }

                if wavelength_cm == T::zero() {
                    return Err(SpectralTransformError::Singular);
                }

                let value = tau / wavelength_cm;
                let first = -tau * length_scale / wavelength_cm.powi(2);
                let second =
                    (T::one() + T::one()) * tau * length_scale.powi(2) / wavelength_cm.powi(3);

                Ok(TransformJet::new(value, first, second))
            }
        }
    }

    #[allow(clippy::wrong_self_convention)]
    pub(crate) fn from_canonical<T>(
        self,
        vacuum_angular_wavenumber: T,
    ) -> Result<TransformJet<T>, SpectralTransformError>
    where
        T: Float + FromPrimitive,
    {
        if !vacuum_angular_wavenumber.is_finite() {
            return Err(SpectralTransformError::NonFinite);
        }

        let tau: T = conversion_constant(std::f64::consts::TAU);
        let c: T = conversion_constant(SPEED_OF_LIGHT_CM_PER_S);

        match self {
            Self::VacuumAngularWavenumber(unit) => {
                let inverse_scale = T::one() / unit.scale();

                Ok(TransformJet::new(
                    vacuum_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::VacuumWavenumber(unit) => {
                let inverse_scale = T::one() / (tau * unit.scale());

                Ok(TransformJet::new(
                    vacuum_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::Frequency(unit) => {
                let inverse_scale = c / (tau * unit.scale());

                Ok(TransformJet::new(
                    vacuum_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::AngularFrequency(unit) => {
                let inverse_scale = c / unit.scale();

                Ok(TransformJet::new(
                    vacuum_angular_wavenumber * inverse_scale,
                    inverse_scale,
                    T::zero(),
                ))
            }

            Self::VacuumWavelength(unit) => {
                if !vacuum_angular_wavenumber.is_finite() {
                    return Err(SpectralTransformError::NonFinite);
                }

                if vacuum_angular_wavenumber == T::zero() {
                    return Err(SpectralTransformError::Singular);
                }

                let length_scale: T = unit.scale();

                let value = tau / (length_scale * vacuum_angular_wavenumber);

                let first = -tau / (length_scale * vacuum_angular_wavenumber.powi(2));

                let second = (T::one() + T::one()) * tau
                    / (length_scale * vacuum_angular_wavenumber.powi(3));

                Ok(TransformJet::new(value, first, second))
            }
        }
    }

    pub(crate) fn metadata(&self) -> AxisMetadata {
        AxisMetadata {
            quantity: match self {
                Self::VacuumAngularWavenumber(_) => AxisQuantity::VacuumAngularWavenumber,
                Self::VacuumWavenumber(_) => AxisQuantity::VacuumWavenumber,
                Self::Frequency(_) => AxisQuantity::Frequency,
                Self::AngularFrequency(_) => AxisQuantity::AngularFrequency,
                Self::VacuumWavelength(_) => AxisQuantity::VacuumWavelength,
            },
            quantity_name: self.quantity_name(),
            quantity_symbol: self.quantity_symbol(),
            unit_symbol: Some(self.unit_symbol()),
            label: self.axis_label(),
        }
    }
}

impl CoordinateLabel for SpectralCoordinate {
    fn quantity_name(&self) -> &'static str {
        match self {
            Self::VacuumAngularWavenumber(_) => "Vacuum angular wavenumber",
            Self::VacuumWavenumber(_) => "Vacuum wavenumber",
            Self::Frequency(_) => "Frequency",
            Self::AngularFrequency(_) => "Angular frequency",
            Self::VacuumWavelength(_) => "Vacuum wavelength",
        }
    }

    fn quantity_symbol(&self) -> &'static str {
        match self {
            Self::VacuumAngularWavenumber(_) => "k₀",
            Self::VacuumWavenumber(_) => "ν̃",
            Self::Frequency(_) => "f",
            Self::AngularFrequency(_) => "ω",
            Self::VacuumWavelength(_) => "λ₀",
        }
    }

    fn unit_symbol(&self) -> &'static str {
        match self {
            Self::VacuumAngularWavenumber(unit) | Self::VacuumWavenumber(unit) => unit.symbol(),
            Self::Frequency(unit) => unit.symbol(),
            Self::AngularFrequency(unit) => unit.symbol(),
            Self::VacuumWavelength(unit) => unit.symbol(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-13_f64.max(1.0e-12 * expected.abs());

        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:e}, got {actual:e}; \
             absolute error was {:e}",
            (actual - expected).abs(),
        );
    }

    fn all_coordinates() -> [SpectralCoordinate; 9] {
        [
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerMetre),
            SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre),
            SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerMicrometre),
            SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
            SpectralCoordinate::Frequency(FrequencyUnit::Hertz),
            SpectralCoordinate::AngularFrequency(AngularFrequencyUnit::TeraradiansPerSecond),
            SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
            SpectralCoordinate::VacuumWavelength(LengthUnit::Nanometre),
        ]
    }

    #[test]
    fn vacuum_angular_wavenumber_in_inverse_centimetres_is_identity() {
        let coordinate =
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre);

        let jet = coordinate.to_canonical(12.5_f64).unwrap();

        assert_close(jet.value, 12.5);
        assert_close(jet.first, 1.0);
        assert_eq!(jet.second, 0.0);
    }

    #[test]
    fn vacuum_angular_wavenumber_units_are_scaled_to_inverse_centimetres() {
        let coordinate = SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerMetre);

        let jet = coordinate.to_canonical(500.0_f64).unwrap();

        assert_close(jet.value, 5.0);
        assert_close(jet.first, 1.0e-2);
        assert_eq!(jet.second, 0.0);
    }

    #[test]
    fn vacuum_wavenumber_is_multiplied_by_two_pi() {
        let coordinate = SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre);

        let jet = coordinate.to_canonical(1.0_f64).unwrap();

        assert_close(jet.value, std::f64::consts::TAU);
        assert_close(jet.first, std::f64::consts::TAU);
        assert_eq!(jet.second, 0.0);
    }

    #[test]
    fn inverse_micrometre_vacuum_wavenumber_converts_correctly() {
        let coordinate = SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerMicrometre);

        let jet = coordinate.to_canonical(0.1_f64).unwrap();

        let expected = std::f64::consts::TAU * 0.1 * 1.0e4;

        assert_close(jet.value, expected);
        assert_close(jet.first, std::f64::consts::TAU * 1.0e4);
        assert_eq!(jet.second, 0.0);
    }

    #[test]
    fn frequency_converts_to_vacuum_angular_wavenumber() {
        let coordinate = SpectralCoordinate::Frequency(FrequencyUnit::Terahertz);

        let jet = coordinate.to_canonical(30.0_f64).unwrap();

        let expected = std::f64::consts::TAU * 30.0e12 / SPEED_OF_LIGHT_CM_PER_S;

        assert_close(jet.value, expected);
        assert_close(
            jet.first,
            std::f64::consts::TAU * 1.0e12 / SPEED_OF_LIGHT_CM_PER_S,
        );
        assert_eq!(jet.second, 0.0);
    }

    #[test]
    fn angular_frequency_converts_to_vacuum_angular_wavenumber() {
        let coordinate =
            SpectralCoordinate::AngularFrequency(AngularFrequencyUnit::TeraradiansPerSecond);

        let jet = coordinate
            .to_canonical(30.0 * std::f64::consts::TAU)
            .unwrap();

        let expected = std::f64::consts::TAU * 30.0e12 / SPEED_OF_LIGHT_CM_PER_S;

        assert_close(jet.value, expected);
    }

    #[test]
    fn wavelength_converts_to_vacuum_angular_wavenumber() {
        let coordinate = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        let jet = coordinate.to_canonical(10.0_f64).unwrap();

        let wavelength_cm = 10.0e-4;
        let expected = std::f64::consts::TAU / wavelength_cm;

        assert_close(jet.value, expected);
        assert_close(
            jet.first,
            -std::f64::consts::TAU * 1.0e-4 / wavelength_cm.powi(2),
        );
        assert_close(
            jet.second,
            2.0 * std::f64::consts::TAU * 1.0e-8 / wavelength_cm.powi(3),
        );
    }

    #[test]
    fn vacuum_wavenumber_and_wavelength_are_consistent() {
        let from_wavelength = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre)
            .to_canonical(10.0_f64)
            .unwrap();

        let from_wavenumber =
            SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre)
                .to_canonical(1_000.0_f64)
                .unwrap();

        assert_close(from_wavelength.value, from_wavenumber.value);
    }

    #[test]
    fn frequency_and_wavelength_are_consistent() {
        let wavelength_cm = 10.0e-4;
        let frequency_hz = SPEED_OF_LIGHT_CM_PER_S / wavelength_cm;

        let from_frequency = SpectralCoordinate::Frequency(FrequencyUnit::Hertz)
            .to_canonical(frequency_hz)
            .unwrap();

        let from_wavelength = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre)
            .to_canonical(10.0_f64)
            .unwrap();

        assert_close(from_frequency.value, from_wavelength.value);
    }

    #[test]
    fn ordinary_and_angular_frequency_are_consistent() {
        let frequency_hz = 30.0e12_f64;
        let angular_frequency = std::f64::consts::TAU * frequency_hz;

        let from_frequency = SpectralCoordinate::Frequency(FrequencyUnit::Hertz)
            .to_canonical(frequency_hz)
            .unwrap();

        let from_angular_frequency =
            SpectralCoordinate::AngularFrequency(AngularFrequencyUnit::RadiansPerSecond)
                .to_canonical(angular_frequency)
                .unwrap();

        assert_close(from_frequency.value, from_angular_frequency.value);
    }

    #[test]
    fn all_coordinates_round_trip() {
        for coordinate in all_coordinates() {
            let original = match coordinate {
                SpectralCoordinate::VacuumWavelength(_) => 10.0,
                _ => 12.5,
            };

            let canonical = coordinate.to_canonical(original).unwrap();

            let recovered = coordinate.from_canonical(canonical.value).unwrap();

            assert_close(recovered.value, original);
        }
    }

    #[test]
    fn forward_and_inverse_first_derivatives_are_reciprocal() {
        for coordinate in all_coordinates() {
            let original = match coordinate {
                SpectralCoordinate::VacuumWavelength(_) => 10.0,
                _ => 12.5,
            };

            let forward = coordinate.to_canonical(original).unwrap();

            let inverse = coordinate.from_canonical(forward.value).unwrap();

            assert_close(forward.first * inverse.first, 1.0);
        }
    }

    #[test]
    fn forward_and_inverse_second_derivatives_compose_to_zero() {
        for coordinate in all_coordinates() {
            let original = match coordinate {
                SpectralCoordinate::VacuumWavelength(_) => 10.0,
                _ => 12.5,
            };

            let forward = coordinate.to_canonical(original).unwrap();

            let inverse = coordinate.from_canonical(forward.value).unwrap();

            let composed_second =
                inverse.second * forward.first.powi(2) + inverse.first * forward.second;

            assert_close(composed_second, 0.0);
        }
    }

    #[test]
    fn zero_wavelength_is_rejected() {
        let coordinate = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        let result = coordinate.to_canonical(0.0_f64);

        assert_eq!(result, Err(SpectralTransformError::Singular));
    }

    #[test]
    fn zero_canonical_wavenumber_cannot_be_converted_to_wavelength() {
        let coordinate = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        let result = coordinate.from_canonical(0.0_f64);

        assert_eq!(result, Err(SpectralTransformError::Singular));
    }

    #[test]
    fn non_finite_wavelength_is_rejected() {
        let coordinate = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        assert_eq!(
            coordinate.to_canonical(f64::NAN),
            Err(SpectralTransformError::NonFinite),
        );

        assert_eq!(
            coordinate.to_canonical(f64::INFINITY),
            Err(SpectralTransformError::NonFinite),
        );
    }

    #[test]
    fn coordinate_labels_are_correct() {
        let angular = SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre);

        assert_eq!(angular.quantity_name(), "Vacuum angular wavenumber",);
        assert_eq!(angular.quantity_symbol(), "k₀");
        assert_eq!(angular.unit_symbol(), "cm⁻¹");
        assert_eq!(angular.axis_label(), "k₀ (cm⁻¹)");

        let wavenumber = SpectralCoordinate::VacuumWavenumber(InverseLengthUnit::PerCentimetre);

        assert_eq!(wavenumber.quantity_name(), "Vacuum wavenumber",);
        assert_eq!(wavenumber.quantity_symbol(), "ν̃");
        assert_eq!(wavenumber.axis_label(), "ν̃ (cm⁻¹)");

        let wavelength = SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre);

        assert_eq!(wavelength.quantity_symbol(), "λ₀");
        assert_eq!(wavelength.axis_label(), "λ₀ (µm)");
    }
}
