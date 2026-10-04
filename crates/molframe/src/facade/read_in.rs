//! Reads that run their perception pass on a caller-supplied execution context.

use super::{dispatch_read, enrich_read};
use crate::structure::Structure;
use molframe_core::ExecutionContext;
use molframe_core::diagnostic::{Diagnostic, Findings};
use molframe_core::io::{InputBuffer, ReadOptions};
use std::path::Path;

/// [`read_with_options`](super::read_with_options) on the worker budget,
/// memory budget and cancellation of `context`.
///
/// The context governs bond and secondary-structure perception, the only part
/// of a read that is parallel or large.
///
/// # Errors
///
/// Returns the findings that stopped the read, including a cancellation.
pub fn read_with_options_in(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    let path = path.as_ref();
    let input = InputBuffer::open(path, options.limits).map_err(Findings::from)?;
    let name = path.file_name().and_then(|name| name.to_str());
    read_buffer_in(&input, name, options, context)
}

/// [`read_bytes`](super::read_bytes) on the budgets of `context`.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_bytes_in(
    bytes: Vec<u8>,
    name: Option<&str>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    read_buffer_in(&InputBuffer::from_bytes(bytes), name, options, context)
}

/// [`read_buffer`](super::read_buffer) on the budgets of `context`.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_buffer_in(
    input: &InputBuffer,
    name: Option<&str>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    enrich_read(dispatch_read(input, name, options), options, context)
        .map(|(structure, diagnostics)| (structure.into(), diagnostics))
}
