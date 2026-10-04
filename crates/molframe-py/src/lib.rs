//! Curated Python bindings for the `MolFrame` facade.
//!
//! Python mirrors the stable Rust contract rather than the workspace crate
//! graph. Native storage and kernels remain in Rust; this crate owns only
//! lifetime-safe Python views and small conversion boundaries.

#![deny(unsafe_op_in_unsafe_fn)]

mod analysis_result;
#[cfg(any(feature = "analysis", feature = "validation", feature = "spatial"))]
mod backend;
mod batches;
mod bindings;
mod catalog;
#[cfg(feature = "chemistry")]
mod chemistry;
#[cfg(feature = "compare")]
mod compare;
#[cfg(feature = "crystal")]
mod crystal;
#[cfg(feature = "crystal")]
mod crystal_links;
#[cfg(feature = "crystal")]
mod crystal_reduction;
#[cfg(feature = "crystal")]
mod crystal_statistics;
mod editing;
#[cfg(feature = "analysis")]
mod electrostatics;
mod entities;
mod error;
mod execution;
mod formats;
#[cfg(feature = "analysis")]
mod governed;
#[cfg(feature = "analysis")]
mod governed_structure;
mod hierarchy;
mod interop;
#[cfg(feature = "motif")]
mod motif;
mod native_source;
mod policy;
#[cfg(feature = "query")]
mod query_aliases;
mod query_cache;
mod query_messages;
mod reading;
mod secondary;
#[cfg(feature = "query")]
mod selection_expr;
#[cfg(feature = "sequence")]
mod sequence;
#[cfg(feature = "spatial")]
mod spatial;
mod structure_data;
#[cfg(feature = "surface")]
mod surface;
mod table;
#[cfg(feature = "trajectory")]
mod trajectory;
#[cfg(feature = "validation")]
mod validation;
#[cfg(feature = "analysis")]
mod workflow;

pub use native_source::{NativeAtom, NativeBond, NativeStructureSource, NativeTopology};

#[cfg(test)]
#[path = "surface_tests.rs"]
mod surface_tests;
#[cfg(test)]
#[path = "module_tests.rs"]
mod tests;

mod conversion;
mod module;

pub use conversion::structure_from_python;
