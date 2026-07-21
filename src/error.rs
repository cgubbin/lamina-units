#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
pub enum TransformError {
    #[error("coordinate value must be finite")]
    NonFinite,

    #[error("coordinate value must not be singular")]
    Singular,

    #[error("incident-angle parameterisation requires an incident refractive index")]
    MissingIncidentIndex,

    #[error("incident refractive index must be finite")]
    InvalidIncidentIndex,
}
