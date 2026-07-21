//! Units and coordinate transformations for the TMM ecosystem.
//!
//! `tmm-units` converts user-facing physical coordinates into the canonical
//! coordinates used by `tmm-core`. It also carries the first and second
//! derivatives required to transform sensitivities back into caller units.
//!
//! The canonical backend units are:
//!
//! - length and thickness: centimetres;
//! - vacuum wavenumber: inverse centimetres;
//! - parallel wavenumber: inverse centimetres.

mod angle;
mod error;
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
