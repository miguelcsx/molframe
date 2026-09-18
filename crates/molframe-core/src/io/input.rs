//! Bounded owned, mapped, streamed, and transparently decompressed input.

mod buffer;
mod collect;
mod compression;
mod decode;
mod limits;
mod source;
mod source_file;
mod spill_source;

pub use buffer::{InputBuffer, InputKind};
pub use compression::Compression;
pub use limits::Limits;
pub use source::{ByteWindow, SourceBytes, WindowedFile};
pub use source_file::WindowedSourceFile;
pub use spill_source::SpillWindowedFile;

#[cfg(test)]
use crate::diagnostic::Code;
#[cfg(all(test, feature = "mmap"))]
use buffer::MMAP_MIN_BYTES;
#[cfg(test)]
use collect::check_expansion;
#[cfg(all(test, feature = "mmap"))]
use std::fs::File;

#[cfg(test)]
#[path = "input_tests.rs"]
mod tests;
