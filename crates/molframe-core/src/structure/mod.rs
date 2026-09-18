//! The normalised structure, the views over it, and the invariants it holds.
//!
//! A structure is immutable and shared by reference; editing produces a new one.
//! That is what makes handing out a pointer to the coordinates safe, and it is
//! why a view keeps working after the structure it came from has been superseded
//! — it holds its own reference to the snapshot it was made against.

#[cfg(test)]
mod fixture;

mod data;
mod diff;
mod disorder;
mod edit;
mod extensions;
#[cfg(test)]
mod extensions_tests;
mod handle;
mod handle_helpers;
mod materialize;
mod merge;
mod merge_annotations;
mod overlay;
mod sequence;
mod validate;
mod view;

pub use data::*;
pub use diff::*;
pub use edit::*;
pub use extensions::*;
pub use handle::*;
pub use overlay::*;
pub use sequence::*;
pub use validate::*;
pub use view::*;
