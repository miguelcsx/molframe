//! Chunked columnar storage, the provider machinery, and the snapshot
//! representation the curated [`crate::Structure`] wraps.
//!
//! This is `molframe-core`'s own storage layer, re-exported unchanged. A
//! caller writing a custom format reader or a new provider needs it; a caller
//! reading and analysing structures does not, which is why it lives here and
//! not at the crate root.

pub use molframe_core::chunk::{
    AtomChunk, AtomChunkStats, AtomRecord, ChunkBuilder, ElementMask, Extremes, ParentMapping,
    TARGET_CHUNK_ATOMS,
};
pub use molframe_core::column::{BitVec, EncodedColumn, Presence, ValidityMask};
pub use molframe_core::coords::{CoordinateBlock, CoordinateGeneration};
pub use molframe_core::provider::{
    AtomEndpoint, BondChunk, BondChunkProvider, BondChunkRecord, ChunkDescriptor, ChunkId,
    ChunkLayout, DatasetCatalog, DatasetDescriptor, DatasetId, FrameChunk, FrameChunkProvider,
    LocalRow, LogicalRow, PayloadKind, PropertyChunk, PropertyChunkProvider, PropertyKind,
    PropertyValue, ProviderError, StructureChunk, StructureChunkProvider, TARGET_CHUNK_BONDS,
};
pub use molframe_core::selection::AtomSelection;
pub use molframe_core::structure::{CoordinateStore, Structure, StructureData, StructureView};
pub use molframe_core::symbol::Interner;
