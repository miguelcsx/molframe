//! Arrow and `DLPack` interoperability for molframe structure tables.
//!
//! Raw pointers are confined to audited ABI adapters that retain immutable
//! `Structure` snapshots for the full foreign-buffer lifetime. This boundary
//! is independent from the operating-system mapping boundary in `molframe-mmap`.

#![deny(unsafe_op_in_unsafe_fn)]

mod columnar;
mod dataset;
mod graph;
mod numeric;
mod tensor;

pub use columnar::{
    ArrowStream, AtomTable, BondTable, ChainTable, ExportCost, MolframeExtension, ResidueTable,
    TableFileError, extension_name, write_atom_ipc, write_atom_ipc_with_metadata,
    write_atom_parquet, write_atom_parquet_with_metadata,
};
pub use dataset::{
    Dataset, DatasetError, DatasetFilter, DatasetSplit, DatasetWarning, LoadError, ManifestEntry,
    SplitOptions, SplitRatios, SplitStrategy,
};
pub use graph::{
    EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphError, GraphOptions, MissingFeaturePolicy,
    NodeFeature, NodeLevel, graph,
};
pub use tensor::{DLDataType, DLDevice, DLManagedTensor, DLTensor, DlpackError, DlpackTensor};
