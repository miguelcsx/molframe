//! Audited operating-system memory mapping boundary.

#![deny(unsafe_op_in_unsafe_fn)]

mod mapped;

pub use mapped::MappedFile;
