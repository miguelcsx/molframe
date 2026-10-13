//! Reading a multi-model structure file as a trajectory.
//!
//! A multi-model PDB or mmCIF file is an ensemble of frames over one topology,
//! which is what a `Trajectory` is. This route reads the structure and hands its
//! models over as frames, honouring the same `frames` and `atoms` selections the
//! trajectory containers take. A ragged ensemble has no shared atom axis, so it
//! is refused with the error that names `ragged_models`.

use crate::bindings::PyStructure;
use crate::trajectory::PyTrajectory;
use crate::trajectory_stream::resolve_slice;
use molframe::trajectory::DenseChunk;
use pyo3::prelude::*;
use pyo3::types::PySlice;
use std::path::Path;

/// The structure format a path or explicit name selects, if it is one.
pub(crate) fn structure_format(path: &Path, format: Option<&str>) -> Option<&'static str> {
    let name = if let Some(name) = format {
        name.to_ascii_lowercase()
    } else {
        let mut name = path.file_name()?.to_str()?.to_ascii_lowercase();
        for suffix in [".gz", ".zst"] {
            if let Some(stripped) = name.strip_suffix(suffix) {
                name = stripped.to_owned();
            }
        }
        name.rsplit_once('.')?.1.to_owned()
    };
    match name.as_str() {
        "pdb" | "ent" => Some("pdb"),
        "cif" | "mmcif" => Some("mmcif"),
        "bcif" => Some("bcif"),
        _ => None,
    }
}

/// Reads `path` as a structure and returns its models as trajectory frames.
pub(crate) fn read_models(
    py: Python<'_>,
    path: &Path,
    format: &'static str,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<&[usize]>,
) -> PyResult<PyTrajectory> {
    let structure = py
        .detach(|| molframe::read(path))
        .map_err(|findings| crate::bindings::findings_error(&findings))?;
    let structure = PyStructure::new(structure);
    structure.require_resolved()?;
    let total = structure.inner.model_count();
    let slice = resolve_slice(frames, || Ok(total))?;
    let source_atoms = structure.inner.coordinates().len();
    if let Some(index) = atoms
        .iter()
        .copied()
        .flatten()
        .find(|index| **index >= source_atoms)
    {
        return Err(crate::error::index(format!(
            "atom index {index} is out of range for {source_atoms} atoms"
        )));
    }
    let block = py
        .detach(|| structure.inner.ensemble_coordinates(None))
        .map_err(|finding| crate::error::from_diagnostic(&finding))?;
    let mut chunk = DenseChunk {
        atoms: atoms.map_or(source_atoms, <[usize]>::len),
        ..DenseChunk::default()
    };
    let stop = slice.stop.min(total);
    for index in (slice.start..stop).step_by(slice.step) {
        let Some(positions) = block.model(index) else {
            break;
        };
        let mut frame = molframe::trajectory::Timestep {
            frame: index,
            ..molframe::trajectory::Timestep::default()
        };
        match atoms {
            Some(atoms) => frame
                .positions
                .extend(atoms.iter().filter_map(|atom| positions.get(*atom))),
            None => frame.positions.extend_from_slice(positions),
        }
        chunk.push(&frame).map_err(crate::error::kernel)?;
    }
    PyTrajectory::from_chunk(py, format.to_owned(), chunk)
}
