use crate::{CompactColumn, ModelCategory};

/// A computational model declared by `ModelCIF`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModelDescription<'a> {
    /// Current `ModelCIF` ordinal identifier.
    pub id: &'a str,
    /// Assembly used by the model.
    pub assembly_id: Option<&'a str>,
    /// Human-readable model name.
    pub name: Option<&'a str>,
    /// Controlled model type.
    pub model_type: Option<&'a str>,
}

/// One modeled target entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Target<'a> {
    /// Core entity identifier.
    pub entity_id: &'a str,
    /// `ModelCIF` data identifier.
    pub data_id: Option<&'a str>,
    /// Target origin.
    pub origin: Option<&'a str>,
}

/// One template used during modeling.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Template<'a> {
    /// Template identifier.
    pub id: &'a str,
    /// Target chain identifier.
    pub target_chain_id: Option<&'a str>,
    /// Human-readable name.
    pub name: Option<&'a str>,
    /// Template origin.
    pub origin: Option<&'a str>,
    /// Template entity type.
    pub entity_type: Option<&'a str>,
}

/// One deterministic step in a modeling protocol.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProtocolStep<'a> {
    /// Protocol identifier.
    pub protocol_id: &'a str,
    /// Step identifier within the protocol.
    pub step_id: &'a str,
    /// Controlled method type.
    pub method_type: &'a str,
    /// Optional name.
    pub name: Option<&'a str>,
    /// Optional details.
    pub details: Option<&'a str>,
    /// Software group used by this step.
    pub software_group_id: Option<&'a str>,
}

/// A software entry associated with a `ModelCIF` software group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SoftwareGroup<'a> {
    /// Group identifier.
    pub group_id: &'a str,
    /// Identifier pointing into the core `software` category.
    pub software_id: &'a str,
    /// Optional parameter group.
    pub parameter_group_id: Option<&'a str>,
}

pub(crate) fn model_description(
    category: &ModelCategory,
    row: usize,
) -> Option<ModelDescription<'_>> {
    Some(ModelDescription {
        id: required_id(category, "ordinal_id", row)?,
        assembly_id: optional_id(category, "assembly_id", row),
        name: optional_id(category, "model_name", row),
        model_type: optional_id(category, "model_type", row),
    })
}

pub(crate) fn target(category: &ModelCategory, row: usize) -> Option<Target<'_>> {
    Some(Target {
        entity_id: required_id(category, "entity_id", row)?,
        data_id: optional_id(category, "data_id", row),
        origin: optional_id(category, "origin", row),
    })
}

pub(crate) fn template(category: &ModelCategory, row: usize) -> Option<Template<'_>> {
    Some(Template {
        id: required_id(category, "template_id", row)?,
        target_chain_id: optional_id(category, "target_asym_id", row),
        name: optional_id(category, "template_name", row),
        origin: optional_id(category, "template_origin", row),
        entity_type: optional_id(category, "template_entity_type", row),
    })
}

pub(crate) fn protocol(category: &ModelCategory, row: usize) -> Option<ProtocolStep<'_>> {
    Some(ProtocolStep {
        protocol_id: required_id(category, "protocol_id", row)?,
        step_id: required_id(category, "step_id", row)?,
        method_type: required_id(category, "method_type", row)?,
        name: optional_id(category, "step_name", row),
        details: optional_id(category, "details", row),
        software_group_id: optional_id(category, "software_group_id", row),
    })
}

pub(crate) fn software(category: &ModelCategory, row: usize) -> Option<SoftwareGroup<'_>> {
    Some(SoftwareGroup {
        group_id: required_id(category, "group_id", row)?,
        software_id: required_id(category, "software_id", row)?,
        parameter_group_id: optional_id(category, "parameter_group_id", row),
    })
}

fn required_id<'a>(category: &'a ModelCategory, item: &str, row: usize) -> Option<&'a str> {
    category.column(item)?.identifier(row)
}

fn optional_id<'a>(category: &'a ModelCategory, item: &str, row: usize) -> Option<&'a str> {
    category
        .column(item)
        .and_then(|column: &CompactColumn| column.identifier(row))
}
