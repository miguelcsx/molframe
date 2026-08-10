use crate::{
    ModelCategory, ModelDescription, ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template,
};
use pdbiox_core::Structure;

/// Stable key for typed `ModelCIF` metadata on a structure.
pub const MODEL_CIF_EXTENSION: &str = "pdbiox.modelcif.v1";

/// Typed and forward-compatible `ModelCIF` contents.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ModelCif {
    /// Every `ma_*` category, including fields newer than this crate.
    pub categories: Vec<ModelCategory>,
    /// Declared computational models.
    pub models: Vec<ModelDescription>,
    /// Modeled targets.
    pub targets: Vec<Target>,
    /// Modeling templates.
    pub templates: Vec<Template>,
    /// Ordered protocol steps.
    pub protocol_steps: Vec<ProtocolStep>,
    /// Software group membership.
    pub software_groups: Vec<SoftwareGroup>,
    /// Typed confidence and quality metrics.
    pub quality: QualityMetrics,
}

impl ModelCif {
    /// Whether the document contained no `ModelCIF` categories.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.categories.is_empty()
    }

    /// Typed confidence metrics.
    #[must_use]
    pub const fn confidence(&self) -> &QualityMetrics {
        &self.quality
    }
}

/// `ModelCIF` access on structures read through the facade.
pub trait ModelCifExt {
    /// Attached `ModelCIF` metadata, when present.
    fn model_cif(&self) -> Option<&ModelCif>;

    /// Typed confidence metrics, when `ModelCIF` metadata is present.
    fn confidence(&self) -> Option<&QualityMetrics> {
        self.model_cif().map(ModelCif::confidence)
    }
}

impl ModelCifExt for Structure {
    fn model_cif(&self) -> Option<&ModelCif> {
        self.extensions().get(MODEL_CIF_EXTENSION)
    }
}
