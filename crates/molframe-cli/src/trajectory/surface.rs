//! Bounded topology/radius projection and frame-ordered surface output.

use crate::exit::Exit;
use crate::report::{Context, Json, RowWriter};
use molframe_core::{Backpressure, BatchDemand, BatchSource, execution::Retained};
use std::path::Path;

pub(crate) fn sasa(
    path: &Path,
    topology: &Path,
    probe: f32,
    samples: u16,
    radius_set: molframe::chemistry::RadiusSet,
    context: Context,
) -> Exit {
    if !super::streaming::supported(path) {
        eprintln!("bounded SASA streaming requires XTC, DCD or TRR");
        return Exit::Refused;
    }
    let bytes = context.execution.memory_budget().bytes() / 4;
    let mut reader = match super::streaming::open(path, bytes, context.execution) {
        Ok(reader) => reader,
        Err(exit) => return exit,
    };
    let radii = match read_radii(topology, reader.n_atoms(), radius_set, context) {
        Ok(radii) => radii,
        Err(exit) => return exit,
    };
    let mut sink = match RowWriter::new(
        context,
        &["frame", "sasa_angstrom_squared", "time_picosecond"],
    ) {
        Ok(sink) => sink,
        Err(error) => {
            eprintln!("output failed: {error}");
            return Exit::Failure;
        }
    };
    let result = molframe::analysis::sasa_stream(
        &mut *reader,
        &radii,
        probe,
        samples,
        context.execution,
        bytes,
        |frame, time, area| {
            let output = if context.is_json() {
                let mut json = Json::new();
                json.number("frame", frame)
                    .number("sasa_angstrom_squared", area);
                if let Some(time) = time {
                    json.number("time_picosecond", time);
                }
                sink.json_record(&json.finish())
            } else {
                let values = [
                    frame.to_string(),
                    area.to_string(),
                    time.map_or_else(String::new, |time| time.to_string()),
                ];
                sink.row(values.iter().map(String::as_str))
            };
            output.map_err(|error| {
                molframe::trajectory::TrajectoryError::SourceIo {
                    format: "SASA output",
                    kind: error.kind(),
                }
                .into()
            })
        },
    );
    if let Err(error) = result {
        eprintln!("trajectory SASA failed: {error}");
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

// Project radii directly from bounded structural batches. The full topology and
// its coordinates are never materialized merely to determine per-atom radii.
fn read_radii(
    path: &Path,
    atoms: usize,
    set: molframe::chemistry::RadiusSet,
    context: Context,
) -> Result<Retained<Vec<f32>>, Exit> {
    let bytes = atoms
        .checked_mul(size_of::<f32>())
        .ok_or(Exit::Consistency)?;
    let reservation = context.execution.try_reserve(bytes).map_err(|error| {
        eprintln!("radius allocation failed: {error}");
        Exit::Consistency
    })?;
    let mut radii = Vec::with_capacity(atoms);
    let options = molframe::ReadOptions::new()
        .mode(context.mode)
        .missing_element_policy(context.missing_element_policy)
        .ambiguous_residue_boundary_policy(context.residue_boundary_policy);
    let mut source =
        molframe::open_structure_batches(path, &options, context.execution).map_err(|error| {
            eprintln!("topology read failed: {error}");
            Exit::Parse
        })?;
    let origin = path.display().to_string();
    loop {
        let demand = BatchDemand::new(8192, 4 * 1024 * 1024);
        match source
            .next_batch(demand, context.execution)
            .map_err(|error| {
                eprintln!("topology read failed: {error}");
                Exit::Parse
            })? {
            Backpressure::Finished => break,
            Backpressure::Pending => {
                eprintln!("topology batch cannot fit the remaining execution budget");
                return Err(Exit::Consistency);
            }
            Backpressure::Ready(lease) => {
                let batch = lease.batch();
                context.findings(batch.diagnostics(), &origin);
                for &element in batch.elements() {
                    if radii.len() >= atoms {
                        eprintln!(
                            "topology has more atoms than the trajectory; supply a single-model topology"
                        );
                        return Err(Exit::Consistency);
                    }
                    let element = molframe::Element::from_atomic_number(element);
                    let Some(radius) = molframe::chemistry::vdw_radius(element, set) else {
                        eprintln!("radius set has no value for topology atom {}", radii.len());
                        return Err(Exit::Indeterminate);
                    };
                    radii.push(radius);
                }
            }
        }
    }
    if radii.len() != atoms {
        eprintln!("topology and trajectory atom counts differ");
        return Err(Exit::Consistency);
    }
    Retained::new(radii, reservation, bytes).map_err(|error| {
        eprintln!("radius allocation failed: {error}");
        Exit::Consistency
    })
}
