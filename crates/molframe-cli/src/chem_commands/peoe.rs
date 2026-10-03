//! Per-atom partial charges.

use crate::args::PeoeArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::chemistry::{ChargeSource, PeoeOptions};

pub(crate) fn peoe(args: &PeoeArguments, context: Context) -> Exit {
    crate::chemistry::with_ccd(
        args.chemistry.ccd.as_deref(),
        args.chemistry.ccd_version.as_deref(),
        context,
        |ccd, version, context| {
            let structure = match open(&args.input, context) {
                Ok(structure) => structure,
                Err(exit) => return exit,
            };
            let provider = match crate::chemistry::load_ccd(ccd, version, context) {
                Ok(provider) => provider,
                Err(exit) => return exit,
            };
            let structure = match crate::chemistry::annotate(&structure, ccd, version, context) {
                Ok(structure) => structure,
                Err(exit) => return exit,
            };
            let charges = match molframe::chemistry::partial_charges(
                structure.engine(),
                &provider,
                PeoeOptions {
                    passes: args.passes,
                    initial_damping: args.initial_damping,
                    damping_factor: args.damping_factor,
                    minimum_electronegativity_difference: args.minimum_difference,
                    profile: args.profile.into(),
                },
            ) {
                Ok(charges) => charges,
                Err(error) => {
                    eprintln!("partial charges failed: {error:?}");
                    return Exit::Consistency;
                }
            };
            match &charges.source {
                ChargeSource::File => eprintln!("charges: read from the file"),
                ChargeSource::Peoe { dictionary, .. } => {
                    eprintln!("charges: PEOE under CCD {}", dictionary.as_str());
                }
            }
            emit_rows(
                context,
                &["atom", "charge"],
                charges.values.into_iter().enumerate(),
                |(atom, charge)| vec![atom.to_string(), charge.to_string()],
            )
        },
    )
}
