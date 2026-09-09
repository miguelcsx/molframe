//! Typed provider failures.

use super::{ChunkId, DatasetId, LogicalRow};

/// Failure to describe or produce an out-of-core chunk.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ProviderError {
    /// A dataset with rows was configured with no chunks.
    #[error("dataset {dataset} has {rows} rows but no chunks")]
    MissingChunks {
        /// Invalid dataset.
        dataset: DatasetId,
        /// Declared logical rows.
        rows: u64,
    },
    /// A chunk size of zero cannot make progress.
    #[error("dataset {dataset} has a zero-row chunk size")]
    ZeroChunkRows {
        /// Invalid dataset.
        dataset: DatasetId,
    },
    /// Dataset arithmetic exceeded its 64-bit identity space.
    #[error("dataset {dataset} {field} overflows 64-bit provider identity")]
    IdentityOverflow {
        /// Dataset whose layout overflowed.
        dataset: DatasetId,
        /// Arithmetic field that overflowed.
        field: &'static str,
    },
    /// Two datasets use the same stable identity.
    #[error("dataset identity {dataset} is duplicated")]
    DuplicateDataset {
        /// Duplicated identity.
        dataset: DatasetId,
    },
    /// Two dataset chunk ranges overlap.
    #[error("chunk {chunk} belongs to more than one dataset")]
    OverlappingChunkRange {
        /// First overlapping chunk.
        chunk: ChunkId,
    },
    /// A requested dataset does not exist in the catalog.
    #[error("dataset {dataset} is not present in the catalog")]
    UnknownDataset {
        /// Missing identity.
        dataset: DatasetId,
    },
    /// A requested chunk does not belong to the provider.
    #[error("chunk {chunk} does not belong to dataset {dataset}")]
    UnknownChunk {
        /// Provider dataset.
        dataset: DatasetId,
        /// Requested chunk.
        chunk: ChunkId,
    },
    /// A logical row does not fit in the requested dataset.
    #[error("logical row {row} lies outside dataset {dataset}")]
    LogicalRowOutOfRange {
        /// Provider dataset.
        dataset: DatasetId,
        /// Invalid row.
        row: LogicalRow,
    },
    /// A local row does not fit in its chunk.
    #[error("local row {row} lies outside a chunk of {rows} rows")]
    LocalRowOutOfRange {
        /// Invalid local row.
        row: u32,
        /// Chunk row count.
        rows: u32,
    },
    /// The supplied payload has a different row count than its descriptor.
    #[error("chunk payload has {actual} rows but its descriptor declares {declared}")]
    PayloadLengthMismatch {
        /// Descriptor row count.
        declared: u32,
        /// Payload row count.
        actual: u32,
    },
    /// A structure adapter cannot find the requested storage chunk.
    #[error("structure storage chunk {index} is unavailable")]
    StorageChunkUnavailable {
        /// Structure-local chunk ordinal.
        index: usize,
    },
    /// A structure has no requested model coordinate block.
    #[error("structure model {model} has no direct coordinate block")]
    ModelUnavailable {
        /// Structure-local model ordinal.
        model: u32,
    },
    /// A requested property is absent.
    #[error("structure has no atom property named {name}")]
    PropertyUnavailable {
        /// Exact property name.
        name: Box<str>,
    },
    /// A property does not align to the structure's atom rows.
    #[error("property {name} has {actual} rows but the structure has {expected}")]
    PropertyLengthMismatch {
        /// Exact property name.
        name: Box<str>,
        /// Structure atom rows.
        expected: u32,
        /// Property rows.
        actual: u32,
    },
    /// The structure does not carry a complete bond topology.
    #[error("structure bond topology is unavailable")]
    BondTopologyUnavailable,
    /// A bond payload range does not fit the retained immutable bond table.
    #[error("bond rows [{first}, {end}) exceed the retained table of {rows} rows")]
    BondRangeUnavailable {
        /// First requested source bond row.
        first: u64,
        /// Exclusive requested source bond row.
        end: u64,
        /// Rows retained by the source table.
        rows: u64,
    },
}
