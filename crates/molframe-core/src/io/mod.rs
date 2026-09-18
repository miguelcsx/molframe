//! Getting bytes in, recognising what they are, and the contract readers meet.

mod collect_structure;
mod format;
mod input;
mod output;
mod structure_batch;
mod structure_batch_buffer;
mod structure_batch_builder;
mod structure_batch_pool;

pub use collect_structure::collect_structure;
pub use format::{
    AmbiguousResidueBoundaryPolicy, Format, MissingElementPolicy, ParseMode, ReadOptions,
    ReadResult, Reader, Select, SelectAll,
};
pub use input::{
    ByteWindow, Compression, InputBuffer, InputKind, Limits, SourceBytes, SpillWindowedFile,
    WindowedFile, WindowedSourceFile,
};
pub use output::{
    DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES, OutputOptions, OutputSink, TextOutput, write_output,
};
pub use structure_batch::{
    BatchContinuity, ContinuityLevel, StructureAtomRecord, StructureBatch, StructureBatchError,
};
pub use structure_batch_buffer::StructureBatchBuffer;
pub use structure_batch_builder::StructureBatchBuilder;
pub use structure_batch_pool::StructureBatchPool;
