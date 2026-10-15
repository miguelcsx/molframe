//! Typed `ModelCIF` metadata and quality metrics over the shared mmCIF document.

#![forbid(unsafe_code)]

mod dense;
mod lower;
mod model;
mod read;
mod write;

pub use dense::{
    ConfidenceError, DensePae, PaeKind, PaeLayout, PaeMatrix, PaePolicy, PaeStorage,
    PrecisionProbe, ResidueIndex, ResidueRef, read_pae_json,
};
pub use lower::lower;
pub(crate) use model::{
    CompactColumn, DictionaryColumn, DictionaryEntry, IntegerColumn, MemoryBudget, NumberColumn,
    PackedIndices, PackedIntegers, ValueRef, identifier_item,
};
pub use model::{
    GlobalMetric, GlobalMetrics, LocalMetric, LocalMetrics, MODEL_CIF_EXTENSION, MetricDefinition,
    MetricDefinitions, ModelCategory, ModelCif, ModelCifError, ModelCifExt, ModelCifOptions,
    ModelCifReadError, ModelDescription, ModelRow, ModelRows, PairwiseMetric, PairwiseMetrics,
    ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template,
};
pub use read::{ModelCifBatchSource, ModelCifProjection, read_compact, read_compact_with_options};
pub use write::{
    ExportCell, ExportTable, write_canonical, write_canonical_to, write_canonical_with_options,
    write_extension_to,
};
