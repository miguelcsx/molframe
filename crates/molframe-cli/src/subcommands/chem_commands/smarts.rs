//! SMARTS substructure matches in a structure or a dictionary component.

use crate::args::SmartsArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::chemistry::{ComponentProvider as _, SmartsMatch, SmartsPattern};
use std::path::Path;

pub(crate) fn smarts(args: &SmartsArguments, context: Context) -> Exit {
    let pattern = match SmartsPattern::parse(&args.pattern) {
        Ok(pattern) => pattern,
        Err(error) => {
            eprintln!(
                "SMARTS pattern rejected at byte {}: {}",
                error.position, error.message
            );
            return Exit::Usage;
        }
    };
    crate::chemistry::with_ccd(
        args.chemistry.ccd.as_deref(),
        args.chemistry.ccd_version.as_deref(),
        context,
        |ccd, version, context| {
            if args.component {
                component(&pattern, &args.target, ccd, version, context)
            } else {
                structure(&pattern, Path::new(&args.target), ccd, version, context)
            }
        },
    )
}

fn emit(matches: Vec<SmartsMatch>, context: Context, label: impl Fn(usize) -> String) -> Exit {
    let rows = matches.into_iter().enumerate().flat_map(|(index, found)| {
        found
            .atom_indices
            .iter()
            .copied()
            .enumerate()
            .map(|(query, atom)| (index, query, atom))
            .collect::<Vec<_>>()
    });
    emit_rows(
        context,
        &["match", "query_atom", "atom"],
        rows,
        |(index, query, atom)| vec![index.to_string(), query.to_string(), label(atom)],
    )
}

fn structure(
    pattern: &SmartsPattern,
    input: &Path,
    ccd: &Path,
    version: &str,
    context: Context,
) -> Exit {
    let opened = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let annotated = match crate::chemistry::annotate(&opened, ccd, version, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    match pattern.find_structure_matches(annotated.engine()) {
        Ok(matches) => emit(matches, context, |atom| atom.to_string()),
        Err(error) => {
            eprintln!("SMARTS search failed: {error}");
            Exit::Consistency
        }
    }
}

fn component(
    pattern: &SmartsPattern,
    id: &str,
    ccd: &Path,
    version: &str,
    context: Context,
) -> Exit {
    let provider = match crate::chemistry::load_ccd(ccd, version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let component = match provider.get(id) {
        Ok(Some(component)) => component,
        Ok(None) => {
            eprintln!("component {id:?} is absent from {}", ccd.display());
            return Exit::Indeterminate;
        }
        Err(finding) => {
            context.findings(&[finding], &ccd.display().to_string());
            return Exit::Consistency;
        }
    };
    let matches = pattern.find_matches(&component);
    emit(matches, context, |atom| match component.atoms.get(atom) {
        Some(found) => found.name.to_string(),
        None => atom.to_string(),
    })
}
