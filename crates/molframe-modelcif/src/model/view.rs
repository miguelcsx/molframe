use crate::model::metadata;
use crate::{
    ModelCategory, ModelDescription, ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template,
};
use molframe_core::Structure;

/// Stable key for typed `ModelCIF` metadata on a structure.
pub const MODEL_CIF_EXTENSION: &str = "molframe.modelcif.v1";

/// Typed and forward-compatible `ModelCIF` contents.
///
/// Every supported and unknown `ma_*` item lives exactly once in compact
/// columnar storage. Typed accessors return borrowed row views rather than
/// cloning identifiers into a second representation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ModelCif {
    /// Every `ma_*` category, including fields newer than this crate.
    pub categories: Vec<ModelCategory>,
}

impl ModelCif {
    /// Whether the document contained no `ModelCIF` categories.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.categories.is_empty()
    }

    /// Looks up one compact category by its name without a leading underscore.
    #[must_use]
    pub fn category(&self, name: &str) -> Option<&ModelCategory> {
        self.categories
            .iter()
            .find(|category| category.name() == name)
    }

    /// Declared computational models in source-row order.
    pub fn models(&self) -> impl Iterator<Item = ModelDescription<'_>> {
        let category = self.category("ma_model_list");
        (0..category.map_or(0, ModelCategory::row_count)).filter_map(move |row| {
            category.and_then(|value| metadata::model_description(value, row))
        })
    }

    /// Modeled targets in source-row order.
    pub fn targets(&self) -> impl Iterator<Item = Target<'_>> {
        let category = self.category("ma_target_entity");
        (0..category.map_or(0, ModelCategory::row_count))
            .filter_map(move |row| category.and_then(|value| metadata::target(value, row)))
    }

    /// Modeling templates in source-row order.
    pub fn templates(&self) -> impl Iterator<Item = Template<'_>> {
        let category = self.category("ma_template_details");
        (0..category.map_or(0, ModelCategory::row_count))
            .filter_map(move |row| category.and_then(|value| metadata::template(value, row)))
    }

    /// Ordered modeling protocol steps.
    pub fn protocol_steps(&self) -> impl Iterator<Item = ProtocolStep<'_>> {
        let category = self.category("ma_protocol_step");
        (0..category.map_or(0, ModelCategory::row_count))
            .filter_map(move |row| category.and_then(|value| metadata::protocol(value, row)))
    }

    /// Software group membership in source-row order.
    pub fn software_groups(&self) -> impl Iterator<Item = SoftwareGroup<'_>> {
        let category = self.category("ma_software_group");
        (0..category.map_or(0, ModelCategory::row_count))
            .filter_map(move |row| category.and_then(|value| metadata::software(value, row)))
    }

    /// Typed confidence metrics.
    #[must_use]
    pub const fn confidence(&self) -> QualityMetrics<'_> {
        QualityMetrics::new(self)
    }
}

/// `ModelCIF` access on structures read through the facade.
pub trait ModelCifExt {
    /// Attached `ModelCIF` metadata, when present.
    fn model_cif(&self) -> Option<&ModelCif>;

    /// Typed confidence metrics, when `ModelCIF` metadata is present.
    fn confidence(&self) -> Option<QualityMetrics<'_>> {
        self.model_cif().map(ModelCif::confidence)
    }
}

impl ModelCifExt for Structure {
    fn model_cif(&self) -> Option<&ModelCif> {
        self.extensions().get(MODEL_CIF_EXTENSION)
    }
}
