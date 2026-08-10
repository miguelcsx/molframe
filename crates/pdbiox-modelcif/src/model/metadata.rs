/// A computational model declared by `ModelCIF`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModelDescription {
    /// Current `ModelCIF` ordinal identifier.
    pub id: Box<str>,
    /// Assembly used by the model.
    pub assembly_id: Option<Box<str>>,
    /// Human-readable model name.
    pub name: Option<Box<str>>,
    /// Controlled model type.
    pub model_type: Option<Box<str>>,
}

/// One modeled target entity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Target {
    /// Core entity identifier.
    pub entity_id: Box<str>,
    /// `ModelCIF` data identifier.
    pub data_id: Option<Box<str>>,
    /// Target origin.
    pub origin: Option<Box<str>>,
}

/// One template used during modeling.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Template {
    /// Template identifier.
    pub id: Box<str>,
    /// Target chain identifier.
    pub target_chain_id: Option<Box<str>>,
    /// Human-readable name.
    pub name: Option<Box<str>>,
    /// Template origin.
    pub origin: Option<Box<str>>,
    /// Template entity type.
    pub entity_type: Option<Box<str>>,
}

/// One deterministic step in a modeling protocol.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProtocolStep {
    /// Protocol identifier.
    pub protocol_id: Box<str>,
    /// Step identifier within the protocol.
    pub step_id: Box<str>,
    /// Controlled method type.
    pub method_type: Box<str>,
    /// Optional name.
    pub name: Option<Box<str>>,
    /// Optional details.
    pub details: Option<Box<str>>,
    /// Software group used by this step.
    pub software_group_id: Option<Box<str>>,
}

/// A software entry associated with a `ModelCIF` software group.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SoftwareGroup {
    /// Group identifier.
    pub group_id: Box<str>,
    /// Identifier pointing into the core `software` category.
    pub software_id: Box<str>,
    /// Optional parameter group.
    pub parameter_group_id: Option<Box<str>>,
}
