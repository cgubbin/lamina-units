#[derive(Clone, Debug, PartialEq, Eq, Copy, Hash)]
pub enum AxisQuantity {
    VacuumAngularWavenumber,
    VacuumWavenumber,
    Frequency,
    AngularFrequency,
    VacuumWavelength,
    ParallelAngularWavenumber,
    ParallelWavenumber,
    EffectiveIndex,
    IncidentAngle,
}

/// Metadata describing one caller-facing coordinate axis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AxisMetadata {
    pub quantity: AxisQuantity,

    /// Human-readable quantity name.
    pub quantity_name: &'static str,

    /// Mathematical quantity symbol.
    pub quantity_symbol: &'static str,

    /// Unit symbol, or `None` for dimensionless coordinates.
    pub unit_symbol: Option<&'static str>,

    /// Complete display label suitable for a plot axis.
    pub label: String,
}
