//! Residue contact map rows.

use crate::args::ContactMapArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::spatial::SpatialBackend;

pub(crate) fn contact_map(args: &ContactMapArguments, context: Context) -> Exit {
    let structure = match open(&args.input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let map = match molframe::analysis::residue_contact_map(
        structure.engine(),
        args.cutoff,
        args.min_separation,
        SpatialBackend::Auto,
        context.execution,
    ) {
        Ok(map) => map,
        Err(error) => {
            eprintln!("contact map failed: {error}");
            return Exit::Usage;
        }
    };
    emit_rows(
        context,
        &["first_residue", "second_residue", "min_distance_angstrom"],
        map.contacts().iter(),
        |contact| {
            vec![
                contact.first.get().to_string(),
                contact.second.get().to_string(),
                contact.min_distance.to_string(),
            ]
        },
    )
}
