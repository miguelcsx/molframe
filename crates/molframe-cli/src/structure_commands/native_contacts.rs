//! Native contact fraction (Q) between two structures.

use crate::args::{AtomMatch, NativeContactArguments};
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::spatial::SpatialBackend;

pub(crate) fn native_contacts(args: &NativeContactArguments, context: Context) -> Exit {
    let reference = match open(&args.reference, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let target = match open(&args.target, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let scored = match args.atom_match {
        AtomMatch::Index => molframe::analysis::native_contact_fraction(
            reference.engine(),
            target.engine(),
            args.cutoff,
            args.contact_tolerance,
            SpatialBackend::Auto,
            context.execution,
        ),
        AtomMatch::Names => molframe::analysis::native_contact_fraction_matched(
            reference.engine(),
            target.engine().positions(),
            &molframe::analysis::match_atoms_by_name(reference.engine(), target.engine()),
            args.cutoff,
            args.contact_tolerance,
            SpatialBackend::Auto,
            context.execution,
        ),
    };
    let result = match scored {
        Ok(result) => result,
        Err(error) => {
            eprintln!("native contacts failed: {error}");
            return Exit::Consistency;
        }
    };
    emit_rows(
        context,
        &["native", "kept", "fraction"],
        [result],
        |contacts| {
            vec![
                contacts.native.to_string(),
                contacts.kept.to_string(),
                contacts.fraction.to_string(),
            ]
        },
    )
}
