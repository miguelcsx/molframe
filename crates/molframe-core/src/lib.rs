//! Core types, storage and analysis contracts for molframe.
//!
//! This crate holds the vocabulary every other molframe crate speaks: the typed
//! positions into a structure's tables, the interned identifiers, the chunked
//! columnar store the coordinates and annotations live in, the immutable
//! [`Structure`] built over it, and the contract types
//! that make an analysis state what it assumed.
//!
//! It has almost no dependencies and no input or output of its own, so a project
//! that wants only the types pays for only the types. Memory mapping and
//! decompression are behind features for the same reason.
//!
//! # Layout
//!
//! - [`index`] — typed positions; a residue position is not an atom position.
//! - [`element`] — element identity, and recovering one from an atom name.
//! - [`span`] — where in a source file something was.
//! - [`diagnostic`] — findings about data, returned alongside success.
//! - [`symbol`] — interned identifiers, so no kernel compares strings.

#![forbid(unsafe_code)]

pub mod annotation;
pub mod bond;
pub mod chunk;
pub mod column;
pub mod contract;
pub mod coords;
pub mod diagnostic;
pub mod element;
pub mod execution;
pub mod hashing;
pub mod index;
pub mod io;
pub mod limits;
pub mod optional;
pub mod parallel;
pub mod provider;
pub mod selection;
pub mod span;
pub mod structure;
pub mod symbol;
pub mod topology;

pub use annotation::{
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn,
    AtomAnnotation, AtomAnnotations, COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION,
    HBOND_ACCEPTOR_ANNOTATION, HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION,
    PLDDT_ANNOTATION, POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
pub use bond::{BondAdjacency, BondOrder, BondProvenance, BondRecord, BondTable, BondTableBuilder};
pub use chunk::{
    AtomChunk, AtomChunkStats, ChunkBuilder, ElementMask, Extremes, ParentMapping,
    TARGET_CHUNK_ATOMS,
};
pub use column::{BitVec, EncodedColumn, Presence, ValidityMask};
pub use contract::{Analysis, AnalysisPolicy, Coverage, Provenance, Status};
pub use coords::{Aabb, CoordinateBlock, CoordinateGeneration};
pub use diagnostic::{
    Class, Code, ContextItem, Diagnostic, Diagnostics, Findings, Kind, Rendered, Severity,
    Strictness,
};
pub use element::Element;
pub use execution::{
    Backpressure, Batch, BatchDemand, BatchLease, BatchSource, CancellationToken, Cancelled,
    ContextError, ExecutionContext, ExecutionContextBuilder, MemoryBudget, MemoryBudgetError,
    MemoryReservation, ScratchPolicy, SpillArtifact, SpillError, SpillFile, SpillReader,
    TempStoragePolicy,
};
pub use index::{
    AtomIndex, BondIndex, ChainIndex, EntityIndex, InstanceId, ModelIndex, ResidueIndex,
};
pub use io::{
    AmbiguousResidueBoundaryPolicy, BatchContinuity, ByteWindow, ContinuityLevel, Format,
    InputBuffer, InputKind, MissingElementPolicy, ParseMode, ReadOptions, ReadResult, SourceBytes,
    SpillWindowedFile, StructureAtomRecord, StructureBatch, StructureBatchBuffer,
    StructureBatchBuilder, StructureBatchError, StructureBatchPool, WindowedFile,
    WindowedSourceFile, collect_structure, write_output,
};
pub use limits::{CapacityError, TableError};
pub use optional::{OptionalI32, OptionalSymbol};
pub use provider::{
    AtomEndpoint, BondChunk, BondChunkProvider, BondChunkRecord, ChunkDescriptor, ChunkId,
    ChunkLayout, DatasetCatalog, DatasetDescriptor, DatasetId, FrameChunk, FrameChunkProvider,
    LocalRow, LogicalRow, PayloadKind, PropertyChunk, PropertyChunkProvider, PropertyKind,
    PropertyValue, ProviderError, StructureChunk, StructureChunkProvider, TARGET_CHUNK_BONDS,
};
pub use selection::AtomSelection;
pub use span::{ByteSpan, Position};
pub use structure::{
    CountDifference, DifferenceError, ExtensionStore, MetadataDifference, Structure, StructureData,
    StructureDifference, StructureDifferenceOptions, StructureView, ValueDifference,
    structure_difference,
};
pub use symbol::{AltId, DictionaryFull, Interner, SymbolId};
pub use topology::{EntityKind, PolymerKind, Topology};
