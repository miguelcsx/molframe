//! DMS command execution over the native system model.

use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use molframe::traj::{DmsSystem, read_dms, write_dms};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) fn info(input: &Path, context: Context) -> Exit {
    let system = match read(input) {
        Ok(system) => system,
        Err(exit) => return exit,
    };
    emit(context, &system);
    Exit::Success
}

pub(crate) fn copy(input: &Path, output: &Path, context: Context) -> Exit {
    let system = match read(input) {
        Ok(system) => system,
        Err(exit) => return exit,
    };
    match write_dms(output, &system) {
        Ok(()) => {
            emit_copy(context, input, output);
            Exit::Success
        }
        Err(error) => {
            eprintln!("could not write {}: {error}", output.display());
            Exit::Refused
        }
    }
}

fn emit_copy(context: Context, input: &Path, output: &Path) {
    if context.is_json() {
        let mut json = Json::new();
        json.text("input", &input.display().to_string())
            .text("output", &output.display().to_string())
            .text("format", "dms");
        context.result(&json.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["input", "output", "format"]);
        let values = [
            input.display().to_string(),
            output.display().to_string(),
            "dms".to_owned(),
        ];
        table.row(values.iter().map(String::as_str));
        context.result(&table.finish());
    } else {
        context.result(&format!(
            "wrote validated DMS system to {}",
            output.display()
        ));
    }
}

fn read(input: &Path) -> Result<DmsSystem, Exit> {
    read_dms(input).map_err(|error| {
        eprintln!("could not read {}: {error}", input.display());
        Exit::Parse
    })
}

fn emit(context: Context, system: &DmsSystem) {
    let particles = system.topology.particles.len();
    let bonds = system.topology.bonds.len();
    let velocities = system.frame.velocities.iter().flatten().count();
    let cell = usize::from(system.frame.cell.is_some());
    if context.is_json() {
        let mut json = Json::new();
        json.text("format", "dms")
            .number("particles", particles)
            .number("bonds", bonds)
            .number("particles_with_velocity", velocities)
            .number("periodic_cell", cell);
        if let Some(version) = system.version {
            json.text("version", &format!("{}.{}", version.major, version.minor));
        }
        context.result(&json.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(
            delimiter,
            &[
                "format",
                "particles",
                "bonds",
                "particles_with_velocity",
                "periodic_cell",
            ],
        );
        let values = [
            "dms".to_owned(),
            particles.to_string(),
            bonds.to_string(),
            velocities.to_string(),
            cell.to_string(),
        ];
        table.row(values.iter().map(String::as_str));
        context.result(&table.finish());
    } else {
        let mut text = String::new();
        let _ = writeln!(text, "format      dms");
        if let Some(version) = system.version {
            let _ = writeln!(text, "version     {}.{}", version.major, version.minor);
        }
        let _ = writeln!(text, "particles   {particles}");
        let _ = writeln!(text, "bonds       {bonds}");
        let _ = writeln!(text, "velocities  {velocities}/{particles} particles");
        let _ = write!(
            text,
            "cell        {}",
            if cell == 1 { "periodic" } else { "absent" }
        );
        context.result(&text);
    }
}
