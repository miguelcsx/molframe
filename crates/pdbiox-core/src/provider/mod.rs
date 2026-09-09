//! Out-of-core provider contracts for renderer and analysis consumers.
//!
//! Catalog metadata scales with the number of logical datasets, while every
//! payload owns only one shared chunk. Parsing and physical I/O remain with the
//! source that implements the provider; consumers request typed chunks by
//! stable global identity.

mod bond_payload;
mod catalog;
mod error;
mod ids;
mod payload;
mod source;

pub use bond_payload::{AtomEndpoint, BondChunk, BondChunkRecord};
pub use catalog::{ChunkLayout, DatasetCatalog, DatasetDescriptor, PayloadKind};
pub use error::ProviderError;
pub use ids::{ChunkId, DatasetId, LocalRow, LogicalRow};
pub use payload::{
    ChunkDescriptor, FrameChunk, PropertyChunk, PropertyKind, PropertyValue, StructureChunk,
};
pub use source::{
    BondChunkProvider, FrameChunkProvider, PropertyChunkProvider, StructureChunkProvider,
    TARGET_CHUNK_BONDS,
};
