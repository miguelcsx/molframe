//! Typed `ModelCIF` data model and structure extension.

mod category;
mod confidence;
mod error;
mod metadata;
mod packed;
mod schema;
mod view;

pub(crate) use category::{
    CompactColumn, DictionaryColumn, DictionaryEntry, IntegerColumn, NumberColumn, ValueRef,
};
pub use category::{ModelCategory, ModelRow, ModelRows};
pub use confidence::{
    GlobalMetric, GlobalMetrics, LocalMetric, LocalMetrics, MetricDefinition, MetricDefinitions,
    PairwiseMetric, PairwiseMetrics, QualityMetrics,
};
pub(crate) use error::MemoryBudget;
pub use error::{ModelCifError, ModelCifOptions, ModelCifReadError};
pub use metadata::{ModelDescription, ProtocolStep, SoftwareGroup, Target, Template};
pub(crate) use packed::{PackedIndices, PackedIntegers};
pub(crate) use schema::identifier_item;
pub use view::{MODEL_CIF_EXTENSION, ModelCif, ModelCifExt};
