//! Atoms of symmetry mates near a structure.

use crate::args::MatesArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::crystal::{CrystalNeighborOptions, SymmetryExt as _};
use molframe::spatial::SpatialBackend;

pub(super) fn mates(args: &MatesArguments, context: Context) -> Exit {
    let structure = match open(&args.input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let Some(symmetry) = structure.engine().symmetry_set() else {
        eprintln!(
            "{} records no space group this catalogue knows, so it has no symmetry mates",
            args.input.display()
        );
        return Exit::Refused;
    };
    let options = CrystalNeighborOptions {
        backend: SpatialBackend::Auto,
        candidate_limit: args.candidate_limit,
        ..CrystalNeighborOptions::default()
    };
    let neighbours = match molframe::crystal::collect_crystal_neighbors(
        structure.engine(),
        symmetry,
        molframe::ModelIndex::new(0),
        args.cutoff,
        options,
        context.execution,
    ) {
        Ok(neighbours) => neighbours,
        Err(finding) => {
            context.findings(&[finding], &args.input.display().to_string());
            return Exit::Consistency;
        }
    };
    // The search also reports the deposited copy's own atoms; those are the
    // molecule, not its crystal contacts.
    let neighbours = neighbours
        .into_iter()
        .filter(|neighbour| neighbour.is_symmetry_mate(symmetry));
    emit_rows(
        context,
        &[
            "source_atom",
            "image_atom",
            "operation",
            "lattice_a",
            "lattice_b",
            "lattice_c",
            "distance_angstrom",
        ],
        neighbours,
        |mate| {
            vec![
                mate.source_atom.get().to_string(),
                mate.image_atom.get().to_string(),
                mate.operation.to_string(),
                mate.lattice[0].to_string(),
                mate.lattice[1].to_string(),
                mate.lattice[2].to_string(),
                mate.distance_squared.sqrt().to_string(),
            ]
        },
    )
}
