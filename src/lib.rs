//! Coordinate transformations and derivative metadata for planar optical systems.
//!
//! This crate converts caller-facing spectral and in-plane coordinates into the
//! canonical coordinates used by a planar transfer-matrix backend.
//!
//! The canonical coordinates are:
//!
//! - vacuum angular wavenumber `k₀`, in `cm⁻¹`;
//! - parallel angular wavenumber `k∥`, in `cm⁻¹`.
//!
//! Callers may instead parameterise a calculation using quantities that are
//! more natural for experiments, simulations, or mode calculations:
//!
//! - wavelength;
//! - frequency or angular frequency;
//! - vacuum wavenumber or vacuum angular wavenumber;
//! - parallel wavenumber;
//! - effective index;
//! - incidence angle in degrees or radians.
//!
//! The crate evaluates the corresponding coordinate transformation together
//! with its first- and second-order derivatives. These derivatives can be used
//! to transform observable derivatives from backend coordinates into the units
//! and parameterisation selected by the caller.
//!
//! # Coupled plane-wave coordinates
//!
//! Spectral and in-plane coordinates are not always independent.
//!
//! For example, effective index is defined by
//!
//! ```text
//! k∥ = n_eff k₀,
//! ```
//!
//! while an incidence angle in an incident medium with refractive index
//! `n_inc(k₀)` gives
//!
//! ```text
//! k∥ = n_inc(k₀) k₀ sin(θ).
//! ```
//!
//! Consequently, the crate treats a plane-wave parameterisation as the coupled
//! transformation
//!
//! ```text
//! (spectral coordinate, in-plane coordinate) -> (k₀, k∥).
//! ```
//!
//! The resulting [`PlaneWaveCoordinateTransform`] contains the value, Jacobian,
//! and Hessian of both canonical coordinates with respect to the two
//! caller-facing coordinates.
//!
//! This distinction is important when transforming derivatives. A spectral
//! derivative taken at fixed angle or fixed effective index is not the same as
//! a backend derivative taken at fixed `k∥`.
//!
//! # Direct in-plane coordinates
//!
//! Coordinates that do not require material information can be transformed in
//! one step.
//!
//! ```
//! use tmm_units::{
//!     FrequencyUnit,
//!     InPlaneCoordinate,
//!     PlaneWaveCoordinates,
//!     SpectralCoordinate,
//!     InverseLengthUnit,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
//!     InPlaneCoordinate::ParallelAngularWavenumber(
//!         InverseLengthUnit::PerCentimetre,
//!     ),
//! );
//!
//! let transform = coordinates
//!     .to_canonical(30.0_f64, 2_500.0)
//!     .unwrap();
//!
//! let k0 = transform.vacuum_angular_wavenumber;
//! let k_parallel = transform.parallel_angular_wavenumber;
//!
//! assert!(k0.value > 0.0);
//! assert_eq!(k_parallel.value, 2_500.0);
//! assert_eq!(k_parallel.spectral_first, 0.0);
//! assert_eq!(k_parallel.in_plane_first, 1.0);
//! ```
//!
//! Here the in-plane coordinate is already the canonical parallel angular
//! wavenumber, so it is independent of the spectral coordinate.
//!
//! # Effective index
//!
//! Effective index is useful when evaluating guided and leaky modes.
//!
//! ```
//! use tmm_units::{
//!     InPlaneCoordinate,
//!     LengthUnit,
//!     PlaneWaveCoordinates,
//!     SpectralCoordinate,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
//!     InPlaneCoordinate::EffectiveIndex,
//! );
//!
//! let transform = coordinates
//!     .to_canonical(10.0_f64, 1.5)
//!     .unwrap();
//!
//! let k0 = transform.vacuum_angular_wavenumber;
//! let k_parallel = transform.parallel_angular_wavenumber;
//!
//! assert!((k_parallel.value - 1.5 * k0.value).abs() < 1.0e-10);
//!
//! // At fixed effective index, changing the spectral coordinate also changes
//! // the parallel wavenumber.
//! assert!(
//!     (k_parallel.spectral_first - 1.5 * k0.spectral_first).abs()
//!         < 1.0e-10
//! );
//!
//! // The derivative with respect to effective index is k₀.
//! assert!(
//!     (k_parallel.in_plane_first - k0.value).abs()
//!         < 1.0e-10
//! );
//! ```
//!
//! # Incidence angle
//!
//! An incidence-angle parameterisation requires the refractive index of the
//! incident medium at the canonical spectral coordinate.
//!
//! The transformation is therefore performed in two stages:
//!
//! 1. transform the caller's spectral coordinate into `k₀`;
//! 2. evaluate the incident material at that `k₀`;
//! 3. resolve the angle into `k∥`.
//!
//! ```
//! use tmm_units::{
//!     AngleUnit,
//!     InPlaneCoordinate,
//!     LengthUnit,
//!     PlaneWaveCoordinates,
//!     RefractiveIndexJet,
//!     SpectralCoordinate,
//!     TransformJet,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
//!     InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
//! );
//!
//! let spectral = coordinates
//!     .transform_spectral(10.0_f64)
//!     .unwrap();
//!
//! // In a real application, these values would be obtained by evaluating the
//! // incident material at `spectral.value`.
//! //
//! // value  = n(k₀)
//! // first  = dn/dk₀
//! // second = d²n/dk₀²
//! let incident_index: RefractiveIndexJet<f64> =
//!     TransformJet::new(1.5, 0.0, 0.0);
//!
//! let transform = coordinates
//!     .to_canonical_with_incident_index(
//!         spectral,
//!         30.0,
//!         incident_index,
//!     )
//!     .unwrap();
//!
//! let k0 = transform.vacuum_angular_wavenumber.value;
//! let k_parallel = transform.parallel_angular_wavenumber.value;
//!
//! let expected = 1.5 * k0 * 0.5;
//! assert!((k_parallel - expected).abs() < 1.0e-10);
//! ```
//!
//! The supplied refractive-index jet must contain derivatives with respect to
//! canonical vacuum angular wavenumber, not with respect to the caller's
//! spectral coordinate.
//!
//! For a dispersive incident medium, the derivatives `dn/dk₀` and `d²n/dk₀²`
//! contribute to the transformed spectral derivatives of `k∥`.
//!
//! # Reusing a spectral transformation
//!
//! Higher-level crates may already need the canonical spectral coordinate for
//! material evaluation. The spectral transformation can be reused to avoid
//! evaluating it more than once.
//!
//! ```
//! use tmm_units::{
//!     InPlaneCoordinate,
//!     LengthUnit,
//!     PlaneWaveCoordinates,
//!     SpectralCoordinate,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
//!     InPlaneCoordinate::EffectiveIndex,
//! );
//!
//! let spectral = coordinates
//!     .transform_spectral(8.0_f64)
//!     .unwrap();
//!
//! let transform = coordinates
//!     .to_canonical_with_spectral(spectral, 2.0)
//!     .unwrap();
//!
//! assert_eq!(
//!     transform.parallel_angular_wavenumber.value,
//!     2.0 * transform.vacuum_angular_wavenumber.value,
//! );
//! ```
//!
//! [`PlaneWaveCoordinates::to_canonical_with_spectral`] supports parameterisations
//! that do not require an incident refractive index. For incidence-angle
//! coordinates, use
//! [`PlaneWaveCoordinates::to_canonical_with_incident_index`] instead.
//!
//! # Coordinate derivatives
//!
//! Each canonical coordinate is represented by a
//! [`BivariateTransformJet`]. If the caller coordinates are `(s, p)`, its
//! fields have the following meaning:
//!
//! | Field | Meaning |
//! | --- | --- |
//! | `value` | `q(s, p)` |
//! | `spectral_first` | `∂q/∂s` |
//! | `in_plane_first` | `∂q/∂p` |
//! | `spectral_second` | `∂²q/∂s²` |
//! | `mixed_second` | `∂²q/(∂s ∂p)` |
//! | `in_plane_second` | `∂²q/∂p²` |
//!
//! These jets provide the Jacobian and coordinate Hessians needed to transform
//! first and second derivatives of backend observables.
//!
//! Suppose a backend observable `F` is a function of the canonical coordinates
//!
//! ```text
//! F = F(k₀, k∥).
//! ```
//!
//! Its derivative with respect to the caller's spectral coordinate is
//!
//! ```text
//! ∂F/∂s
//!     = ∂F/∂k₀ ∂k₀/∂s
//!     + ∂F/∂k∥ ∂k∥/∂s.
//! ```
//!
//! Full second-order conversion additionally requires the canonical mixed
//! derivative
//!
//! ```text
//! ∂²F/(∂k₀ ∂k∥).
//! ```
//!
//! This mixed derivative is generally required for angle and effective-index
//! parameterisations because `k∥` depends on the spectral coordinate.
//!
//! # Axis metadata
//!
//! [`PlaneWaveCoordinates`] also reports metadata for the two caller-facing
//! axes. Plotting crates can use this metadata without inspecting the internal
//! spectral and in-plane coordinate variants.
//!
//! ```
//! use tmm_units::{
//!     AngleUnit,
//!     AxisQuantity,
//!     InPlaneCoordinate,
//!     LengthUnit,
//!     PlaneWaveCoordinates,
//!     SpectralCoordinate,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::VacuumWavelength(LengthUnit::Micrometre),
//!     InPlaneCoordinate::IncidentAngle(AngleUnit::Degree),
//! );
//!
//! let [spectral_axis, in_plane_axis] = coordinates.axes();
//!
//! assert_eq!(
//!     spectral_axis.quantity,
//!     AxisQuantity::VacuumWavelength,
//! );
//! assert_eq!(spectral_axis.label, "λ₀ (µm)");
//!
//! assert_eq!(
//!     in_plane_axis.quantity,
//!     AxisQuantity::IncidentAngle,
//! );
//! assert_eq!(in_plane_axis.label, "θ (°)");
//! ```
//!
//! Dimensionless axes, such as effective index, have no unit symbol:
//!
//! ```
//! use tmm_units::{
//!     AxisQuantity,
//!     FrequencyUnit,
//!     InPlaneCoordinate,
//!     PlaneWaveAxis,
//!     PlaneWaveCoordinates,
//!     SpectralCoordinate,
//! };
//!
//! let coordinates = PlaneWaveCoordinates::new(
//!     SpectralCoordinate::Frequency(FrequencyUnit::Terahertz),
//!     InPlaneCoordinate::EffectiveIndex,
//! );
//!
//! let axis = coordinates.axis_metadata(PlaneWaveAxis::InPlane);
//!
//! assert_eq!(axis.quantity, AxisQuantity::EffectiveIndex);
//! assert_eq!(axis.unit_symbol, None);
//! assert_eq!(axis.label, "n_eff");
//! ```
//!
//! # Units and angular quantities
//!
//! Ordinary and angular quantities use the same underlying dimensional units
//! where appropriate, but remain distinct coordinate variants.
//!
//! For example:
//!
//! ```text
//! vacuum wavenumber          ν̃ = 1 / λ
//! vacuum angular wavenumber  k₀ = 2π / λ
//! ```
//!
//! Both may be expressed in `cm⁻¹`, but their numerical values differ by a
//! factor of `2π`.
//!
//! The same distinction is retained for parallel wavenumber and parallel
//! angular wavenumber.
//!
//! Frequency and angular frequency use separate unit enums because their unit
//! symbols are conventionally distinct:
//!
//! ```text
//! frequency          f  in Hz
//! angular frequency  ω  in rad·s⁻¹
//! ```
//!
//! # Scope
//!
//! This crate performs coordinate algebra and unit conversion. It does not:
//!
//! - evaluate material models;
//! - choose an incident medium;
//! - decide whether a requested angle is physically propagating;
//! - solve the transfer or scattering matrix;
//! - compute optical observables.
//!
//! Higher-level optical crates are responsible for evaluating the incident
//! refractive index and supplying its derivatives when an incidence-angle
//! parameterisation is used.
#![allow(dead_code)]

mod angle;
mod frequency;
mod in_plane;
mod inverse_length;
mod label;
mod length;
mod metadata;
mod plane_wave;
mod private;
mod spectral;
mod transform;

pub use angle::AngleUnit;
pub use frequency::{AngularFrequencyUnit, FrequencyUnit};
pub use in_plane::{InPlaneCoordinate, RefractiveIndexJet};
pub use inverse_length::InverseLengthUnit;
pub use length::LengthUnit;
pub use metadata::{AxisMetadata, AxisQuantity};
pub use plane_wave::{PlaneWaveAxis, PlaneWaveCoordinateTransform, PlaneWaveCoordinates};
pub use spectral::SpectralCoordinate;
pub use transform::TransformJet;
