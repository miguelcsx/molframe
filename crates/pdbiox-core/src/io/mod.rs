//! Getting bytes in, recognising what they are, and the contract readers meet.

mod format;
mod input;
mod output;

pub use format::{
    AmbiguousResidueBoundaryPolicy, Format, MissingElementPolicy, ParseMode, ReadOptions,
    ReadResult, Reader, Select, SelectAll,
};
pub use input::{Compression, InputBuffer, InputKind, Limits};
pub use output::write_output;
