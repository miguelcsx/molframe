//! Incremental contact-count rows from native trajectory readers.

use crate::exit::Exit;
use crate::report::{Context, Json, RowWriter};
use std::path::Path;

pub(crate) fn contacts(path: &Path, cutoff: f32, context: Context) -> Exit {
    if !super::streaming::supported(path) {
        eprintln!("bounded contact streaming requires XTC, DCD or TRR");
        return Exit::Refused;
    }
    // Leave the majority of the shared account for the spatial index and its
    // admitted workers. Reader and frame allowances remain subordinate limits.
    let bytes = context.execution.memory_budget().bytes() / 4;
    let mut reader = match super::streaming::open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    let mut sink = match RowWriter::new(context, &["frame", "contacts", "time_picosecond"]) {
        Ok(sink) => sink,
        Err(error) => {
            eprintln!("output failed: {error}");
            return Exit::Failure;
        }
    };
    let result = molframe::traj::contact_counts_stream(
        &mut *reader,
        cutoff,
        molframe::SpatialSearchOptions::default(),
        context.execution,
        bytes,
        |frame, time, contacts| {
            let output = if context.is_json() {
                let mut json = Json::new();
                json.number("frame", frame).number("contacts", contacts);
                if let Some(time) = time {
                    json.number("time_picosecond", time);
                }
                sink.json_record(&json.finish())
            } else {
                let values = [
                    frame.to_string(),
                    contacts.to_string(),
                    time.map_or_else(String::new, |time| time.to_string()),
                ];
                sink.row(values.iter().map(String::as_str))
            };
            output.map_err(|error| molframe::traj::TrajectoryError::SourceIo {
                format: "contact output",
                kind: error.kind(),
            })
        },
    );
    if let Err(error) = result {
        eprintln!("trajectory contacts failed: {error}");
        return Exit::Consistency;
    }
    match sink.finish() {
        Ok(()) => Exit::Success,
        Err(error) => {
            eprintln!("output failed: {error}");
            Exit::Failure
        }
    }
}
