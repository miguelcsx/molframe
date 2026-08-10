//! Mechanical Python bindings for the Rust facade.
//!
//! Audited pointer use is confined to NumPy/Arrow/DLPack lifetime adapters;
//! scientific kernels remain in Rust. This ABI boundary is separate from the
//! operating-system mapping boundary in `pdbiox-mmap`.

#![deny(unsafe_op_in_unsafe_fn)]

mod analysis;
mod arrow;
mod atom;
mod bonds;
mod chemistry;
mod compatibility;
mod contract;
mod crystallography;
mod dms;
mod edit;
mod errors;
mod external;
mod geometry;
mod graph;
mod hierarchy;
mod index;
mod intrinsic;
mod io;
mod ml;
mod module;
mod query;
mod science;
mod spatial;
mod structure;
mod trajectory;
