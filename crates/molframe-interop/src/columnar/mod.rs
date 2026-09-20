//! Arrow columnar views of molframe structure tables.
//!
//! Atom, residue, chain and bond tables share one shape: bind an immutable
//! [`Structure`][molframe_core::Structure] snapshot, publish a schema with
//! molframe extension and copy-cost metadata, and export either eager record
//! batches or a lazy Arrow C Stream. File writers serialise the atom table
//! chunk at a time.

mod atoms;
mod bond;
mod chain;
mod contact;
mod extension;
mod ipc;
mod owner;
mod parquet;
mod residue;
mod stream;
mod table;
mod table_file;

pub use atoms::AtomTable;
pub use bond::BondTable;
pub use chain::ChainTable;
pub use contact::ContactArrowTable;
pub use extension::{ExportCost, MolframeExtension, extension_name};
pub use ipc::{write_atom_ipc, write_atom_ipc_with_metadata};
pub use parquet::{write_atom_parquet, write_atom_parquet_with_metadata};
pub use residue::ResidueTable;
pub use stream::ArrowStream;
pub use table_file::TableFileError;

#[cfg(test)]
#[path = "table_tests.rs"]
mod tests;
