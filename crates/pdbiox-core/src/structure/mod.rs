//! The normalised structure, the views over it, and the invariants it holds.
//!
//! A structure is immutable and shared by reference; editing produces a new one.
//! That is what makes handing out a pointer to the coordinates safe, and it is
//! why a view keeps working after the structure it came from has been superseded
//! — it holds its own reference to the snapshot it was made against.

#[cfg(test)]
mod fixture;

mod data;
mod handle;
mod validate;
mod view;

pub use data::{CoordinateStore, EntryMetadata, Structure, StructureData, UnitCell};
pub use handle::{AtomRef, ChainRef, ModelRef, ResidueRef};
pub use validate::validate;
pub use view::StructureView;
