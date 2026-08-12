//! Units and physical quantities used by the Lamina ecosystem.
//!
//! This crate provides unit definitions, labels, and conversions shared by
//! Lamina crates.
//!
//! Higher-level optical coordinate compilation belongs in `lamina-core`.

mod constants;
mod label;
mod metadata;
mod private;
mod quantities;
mod units;

pub use label::{CoordinateLabel, UnitLabel};
pub use metadata::{AxisMetadata, AxisQuantity};

pub use quantities::{
    Angle, AngularFrequency, Frequency, Length, SpectralQuantity, SpectralQuantityExt,
    VacuumAngularWavenumber, VacuumAngularWavenumberSquared, VacuumWavenumber, Wavelength,
};
pub use units::{AngleUnit, AngularFrequencyUnit, FrequencyUnit, InverseLengthUnit, LengthUnit};
