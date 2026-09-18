//! Incremental CLI trajectory summaries and RMSD output.

use crate::exit::Exit;
use crate::report::{Context, Json, RowWriter};
use molframe::traj::{
    FrameAlignment, Timestep, TrajectoryFormat, TrajectoryReader, TrajectoryReaderOptions,
    read_trajectory_in, rmsd_stream, run_analysis_stream,
};
use std::path::Path;

pub(super) fn supported(path: &Path) -> bool {
    matches!(
        TrajectoryFormat::infer(path),
        Some(TrajectoryFormat::Xtc | TrajectoryFormat::Dcd | TrajectoryFormat::Trr)
    )
}

pub(super) fn open(
    path: &Path,
    bytes: usize,
    context: &molframe::ExecutionContext,
) -> Result<Box<dyn TrajectoryReader>, Exit> {
    read_trajectory_in(
        path,
        &TrajectoryReaderOptions {
            memory_limit_bytes: bytes,
            ..TrajectoryReaderOptions::default()
        },
        context,
    )
    .map_err(|error| {
        eprintln!("trajectory read failed: {error}");
        Exit::Consistency
    })
}

pub(super) fn info(path: &Path, topology: Option<&Path>, context: Context) -> Exit {
    let bytes = context.execution.memory_budget().bytes() / 2;
    let mut reader = match open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    if let Some(path) = topology {
        let structure = match crate::commands::open(path, context) {
            Ok(value) => value,
            Err(exit) => return exit,
        };
        if structure.atom_count() as usize != reader.n_atoms() {
            eprintln!("topology and trajectory atom counts differ");
            return Exit::Consistency;
        }
    }
    let mut summary = super::summary::Summary {
        format: reader.format(),
        atoms: reader.n_atoms(),
        frames: 0,
        first_time: None,
        last_time: None,
        frames_with_cell: 0,
        frames_with_velocities: 0,
        frames_with_forces: 0,
    };
    let result = run_analysis_stream(&mut *reader, context.execution, bytes, |frame, _| {
        if summary.frames == 0 {
            summary.first_time = frame.time;
        }
        summary.last_time = frame.time;
        summary.frames += 1;
        summary.frames_with_cell += usize::from(frame.cell.is_some());
        summary.frames_with_velocities += usize::from(frame.velocities.is_some());
        summary.frames_with_forces += usize::from(frame.forces.is_some());
        Ok::<(), molframe::traj::TrajectoryError>(())
    });
    if let Err(error) = result {
        eprintln!("trajectory read failed: {error}");
        return Exit::Consistency;
    }
    super::command::emit_summary(context, &summary);
    Exit::Success
}

pub(super) fn rmsd(path: &Path, reference_index: usize, no_fit: bool, context: Context) -> Exit {
    let bytes = context.execution.memory_budget().bytes() / 3;
    let mut reader = match open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    let _reference_reservation = match context.execution.try_reserve(bytes) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("reference allocation failed: {error}");
            return Exit::Consistency;
        }
    };
    let mut reference = Timestep::default();
    for _ in 0..=reference_index {
        match reader.read_next_bounded(&mut reference, bytes) {
            Ok(true) => {}
            Ok(false) => {
                eprintln!("reference frame is absent");
                return Exit::Consistency;
            }
            Err(error) => {
                eprintln!("reference read failed: {error}");
                return Exit::Consistency;
            }
        }
    }
    drop(reader);
    let mut reader = match open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    let mut sink = match RowWriter::new(context, &["frame", "rmsd_angstrom", "time_picosecond"]) {
        Ok(sink) => sink,
        Err(error) => {
            eprintln!("output failed: {error}");
            return Exit::Failure;
        }
    };
    let alignment = if no_fit {
        FrameAlignment::None
    } else {
        FrameAlignment::Rigid
    };
    let result = rmsd_stream(
        &mut *reader,
        &reference.positions,
        alignment,
        context.execution,
        bytes,
        |index, time, value| {
            let output = if context.is_json() {
                let mut json = Json::new();
                json.number("frame", index).number("rmsd_angstrom", value);
                if let Some(time) = time {
                    json.number("time_picosecond", time);
                }
                sink.json_record(&json.finish())
            } else {
                let values = [
                    index.to_string(),
                    value.to_string(),
                    time.map_or_else(String::new, |time| time.to_string()),
                ];
                sink.row(values.iter().map(String::as_str))
            };
            output.map_err(|error| molframe::traj::TrajectoryError::SourceIo {
                format: "result output",
                kind: error.kind(),
            })
        },
    );
    if let Err(error) = result {
        eprintln!("trajectory RMSD failed: {error}");
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
