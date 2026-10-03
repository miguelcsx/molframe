//! Native contact fraction (Q) between two structures.

use crate::args::NativeContactArguments;
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
    let result = match molframe::analysis::native_contact_fraction(
        reference.engine(),
        target.engine(),
        args.cutoff,
        args.retention,
        SpatialBackend::Auto,
        context.execution,
    ) {
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
