//! Trajectory command execution over the native dispatcher.

use super::summary::{Summary, canonical_length_unit, canonical_time_unit, format_name};
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use pdbiox::traj::{
    FrameAlignment, GsdOptions, TrajectoryData, TrajectoryFormat, TrajectoryReadOptions,
    TrajectoryWriteOptions, TrzWriteOptions, read_trajectory_materialized, rmsd_to_reference,
    write_trajectory,
};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) fn info(
    input: &Path,
    topology: Option<&Path>,
    gsd_length_scale: Option<f64>,
    context: Context,
) -> Exit {
    if super::streaming::supported(input) {
        return super::streaming::info(input, topology, context);
    }
    let data = match read(input, gsd_length_scale) {
        Ok(data) => data,
        Err(exit) => return exit,
    };
    if let Some(path) = topology {
        let structure = match open(path, context) {
            Ok(structure) => structure,
            Err(exit) => return exit,
        };
        let atoms = data.frames.first().map_or(0, |frame| frame.positions.len());
        if structure.atom_count() as usize != atoms {
            eprintln!(
                "topology has {} atoms, but trajectory has {atoms}",
                structure.atom_count()
            );
            return Exit::Consistency;
        }
    }
    emit_summary(context, &Summary::from_data(&data));
    Exit::Success
}

pub(crate) fn convert(
    input: &Path,
    output: &Path,
    stride: usize,
    gsd_length_scale: Option<f64>,
    trz_title: Option<&str>,
    context: Context,
) -> Exit {
    let mut data = match read(input, gsd_length_scale) {
        Ok(data) => data,
        Err(exit) => return exit,
    };
    apply_stride(&mut data, stride);
    write_data(input, output, &data, gsd_length_scale, trz_title, context)
}

pub(crate) fn extract(
    input: &Path,
    output: &Path,
    frames: &[usize],
    gsd_length_scale: Option<f64>,
    trz_title: Option<&str>,
    context: Context,
) -> Exit {
    let data = match read(input, gsd_length_scale) {
        Ok(data) => data,
        Err(exit) => return exit,
    };
    let selected = match data.select_frames(frames) {
        Ok(selected) => selected,
        Err(error) => {
            eprintln!("could not select frames: {error}");
            return Exit::Usage;
        }
    };
    write_data(
        input,
        output,
        &selected,
        gsd_length_scale,
        trz_title,
        context,
    )
}

pub(crate) fn rmsd(
    input: &Path,
    reference: usize,
    no_fit: bool,
    gsd_length_scale: Option<f64>,
    context: Context,
) -> Exit {
    if super::streaming::supported(input) {
        return super::streaming::rmsd(input, reference, no_fit, context);
    }
    let data = match read(input, gsd_length_scale) {
        Ok(data) => data,
        Err(exit) => return exit,
    };
    let alignment = if no_fit {
        FrameAlignment::None
    } else {
        FrameAlignment::Rigid
    };
    let values = match rmsd_to_reference(&data.frames, reference, alignment) {
        Ok(values) => values,
        Err(error) => {
            eprintln!("trajectory RMSD failed: {error}");
            return Exit::Consistency;
        }
    };
    emit_rmsd(context, &data, &values);
    Exit::Success
}

fn emit_rmsd(context: Context, data: &TrajectoryData, values: &[f64]) {
    if context.is_json() {
        let objects = data
            .frames
            .iter()
            .zip(values)
            .map(|(frame, value)| {
                let mut json = Json::new();
                json.number("frame", frame.frame)
                    .number("rmsd_angstrom", value);
                if let Some(time) = frame.time {
                    json.number("time_picosecond", time);
                }
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else {
        let delimiter = match context.delimiter() {
            Some(value) => value,
            None => '\t',
        };
        let mut table = Table::new(delimiter, &["frame", "time_picosecond", "rmsd_angstrom"]);
        for (frame, value) in data.frames.iter().zip(values) {
            let values = [
                frame.frame.to_string(),
                optional_number(frame.time),
                value.to_string(),
            ];
            table.row(values.iter().map(String::as_str));
        }
        context.result(&table.finish());
    }
}

fn write_data(
    input: &Path,
    output: &Path,
    data: &TrajectoryData,
    gsd_length_scale: Option<f64>,
    trz_title: Option<&str>,
    context: Context,
) -> Exit {
    let gsd = match gsd_options(gsd_length_scale) {
        Ok(options) => options,
        Err(exit) => return exit,
    };
    let options = TrajectoryWriteOptions {
        gsd,
        trz: trz_title.map(|title| TrzWriteOptions {
            title: title.into(),
        }),
        ..TrajectoryWriteOptions::default()
    };
    if data.format != TrajectoryFormat::Trz
        && TrajectoryFormat::infer(output) == Some(TrajectoryFormat::Trz)
        && options.trz.is_none()
    {
        eprintln!("conversion to TRZ requires --trz-title because the source has no TRZ title");
        return Exit::Usage;
    }
    match write_trajectory(output, data, &options) {
        Ok(()) => {
            let summary = Summary::from_data(data);
            emit_conversion(context, input, output, &summary);
            Exit::Success
        }
        Err(error) => {
            eprintln!("could not write {}: {error}", output.display());
            Exit::Refused
        }
    }
}

fn read(path: &Path, gsd_length_scale: Option<f64>) -> Result<TrajectoryData, Exit> {
    let options = TrajectoryReadOptions {
        gsd: gsd_options(gsd_length_scale)?,
        ..TrajectoryReadOptions::default()
    };
    read_trajectory_materialized(path, &options).map_err(|error| {
        eprintln!("could not read {}: {error}", path.display());
        Exit::Parse
    })
}

fn gsd_options(scale: Option<f64>) -> Result<Option<GsdOptions>, Exit> {
    scale
        .map(|value| {
            GsdOptions::new(value).map_err(|error| {
                eprintln!("invalid --gsd-length-scale: {error}");
                Exit::Usage
            })
        })
        .transpose()
}

pub(super) fn apply_stride(data: &mut TrajectoryData, stride: usize) {
    data.frames = data.frames.iter().step_by(stride).cloned().collect();
    if let Some(steps) = &mut data.metadata.steps {
        *steps = steps.iter().step_by(stride).copied().collect();
    }
}

pub(super) fn emit_summary(context: Context, summary: &Summary) {
    if context.is_json() {
        context.result(&summary_json(summary));
    } else if let Some(delimiter) = context.delimiter() {
        context.result(&summary_table(summary, delimiter));
    } else {
        context.result(&summary_text(summary));
    }
}

fn summary_json(summary: &Summary) -> String {
    let mut json = Json::new();
    json.text("format", summary.format)
        .number("frames", summary.frames)
        .number("atoms", summary.atoms)
        .text("length_unit", canonical_length_unit())
        .text("time_unit", canonical_time_unit())
        .number("frames_with_cell", summary.frames_with_cell)
        .number("frames_with_velocities", summary.frames_with_velocities)
        .number("frames_with_forces", summary.frames_with_forces);
    if let Some(value) = summary.first_time {
        json.number("first_time", value);
    }
    if let Some(value) = summary.last_time {
        json.number("last_time", value);
    }
    json.finish()
}

fn summary_table(summary: &Summary, delimiter: char) -> String {
    let mut table = Table::new(
        delimiter,
        &[
            "format",
            "frames",
            "atoms",
            "first_time_ps",
            "last_time_ps",
            "frames_with_cell",
            "frames_with_velocities",
            "frames_with_forces",
        ],
    );
    let values = [
        summary.format.to_owned(),
        summary.frames.to_string(),
        summary.atoms.to_string(),
        optional_number(summary.first_time),
        optional_number(summary.last_time),
        summary.frames_with_cell.to_string(),
        summary.frames_with_velocities.to_string(),
        summary.frames_with_forces.to_string(),
    ];
    table.row(values.iter().map(String::as_str));
    table.finish()
}

fn summary_text(summary: &Summary) -> String {
    let mut text = String::new();
    let _ = writeln!(text, "format      {}", summary.format);
    let _ = writeln!(text, "frames      {}", summary.frames);
    let _ = writeln!(text, "atoms       {}", summary.atoms);
    let _ = writeln!(text, "coordinates {}", canonical_length_unit());
    if let (Some(first), Some(last)) = (summary.first_time, summary.last_time) {
        let _ = writeln!(
            text,
            "time        {first} to {last} {}",
            canonical_time_unit()
        );
    }
    let _ = writeln!(
        text,
        "cells       {}/{} frames",
        summary.frames_with_cell, summary.frames
    );
    let _ = writeln!(
        text,
        "velocities  {}/{} frames",
        summary.frames_with_velocities, summary.frames
    );
    let _ = write!(
        text,
        "forces      {}/{} frames",
        summary.frames_with_forces, summary.frames
    );
    text
}

fn optional_number(value: Option<f64>) -> String {
    value.map_or_else(String::new, |number| number.to_string())
}

fn emit_conversion(context: Context, input: &Path, output: &Path, summary: &Summary) {
    if context.is_json() {
        let mut json = Json::new();
        json.text("input", &input.display().to_string())
            .text("output", &output.display().to_string())
            .text("format", format_name_from_path(output))
            .number("frames", summary.frames)
            .number("atoms", summary.atoms);
        context.result(&json.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["input", "output", "format", "frames", "atoms"]);
        let values = [
            input.display().to_string(),
            output.display().to_string(),
            format_name_from_path(output).to_owned(),
            summary.frames.to_string(),
            summary.atoms.to_string(),
        ];
        table.row(values.iter().map(String::as_str));
        context.result(&table.finish());
    } else {
        context.result(&format!(
            "wrote {} frames with {} atoms to {} ({})",
            summary.frames,
            summary.atoms,
            output.display(),
            format_name_from_path(output)
        ));
    }
}

fn format_name_from_path(path: &Path) -> &'static str {
    pdbiox::traj::TrajectoryFormat::infer(path).map_or("unknown", format_name)
}
