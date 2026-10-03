//! Per-atom fluctuation rows from a bounded trajectory stream.

use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use std::path::Path;

pub(crate) fn rmsf(path: &Path, context: Context) -> Exit {
    if !super::streaming::supported(path) {
        eprintln!("bounded RMSF streaming requires XTC, DCD or TRR");
        return Exit::Refused;
    }
    let bytes = context.execution.memory_budget().bytes() / 3;
    let mut reader = match super::streaming::open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    let values = match molframe::trajectory::rmsf_stream(&mut *reader, context.execution, bytes) {
        Ok(values) => values,
        Err(error) => {
            eprintln!("trajectory RMSF failed: {error}");
            return Exit::Consistency;
        }
    };
    emit_rows(
        context,
        &["atom", "rmsf_angstrom"],
        values.iter().enumerate(),
        |(atom, value)| vec![atom.to_string(), value.to_string()],
    )
}
