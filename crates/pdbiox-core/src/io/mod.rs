//! Getting bytes in, recognising what they are, and the contract readers meet.

mod format;
mod input;

pub use format::{Format, ParseMode, ReadOptions, ReadResult, Reader, Select, SelectAll};
pub use input::{Compression, InputBuffer, Limits};
