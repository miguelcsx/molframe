//! Format dispatch and the shared reader/writer contract.

mod detect;
mod options;

#[cfg(test)]
use super::input::InputBuffer;
#[cfg(test)]
use crate::diagnostic::{Code, Strictness};

pub use detect::Format;
pub use options::{
    AmbiguousResidueBoundaryPolicy, MissingElementPolicy, ParseMode, ReadOptions, ReadResult,
    Reader, Select, SelectAll,
};

#[cfg(test)]
#[path = "../format_tests.rs"]
mod tests;
