//! Thin rendering over native analysis kernels.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, RowWriter};
use molframe::QueryStructure as _;
use molframe::{AtomIndex, RadiusSet, SpatialBackend};
use std::path::Path;

pub(crate) fn contacts(
    input: &Path,
    between: Option<&[String]>,
    cutoff: f32,
    context: Context,
) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let filters = match contact_filters(&structure, between, input, context) {
        Ok(filters) => filters,
        Err(exit) => return exit,
    };
    let mut output =
        match RowWriter::new(context, &["first_atom", "second_atom", "distance_angstrom"]) {
            Ok(output) => output,
            Err(error) => return output_error(&error),
        };
    let mut write_error = None;
    let mut emit = |contact: molframe::analysis::Contact| {
        if write_error.is_some() {
            return;
        }
        let first = contact.first.get().to_string();
        let second = contact.second.get().to_string();
        let distance = contact.distance.to_string();
        write_error = output
            .row([first.as_str(), second.as_str(), distance.as_str()])
            .err();
    };
    let searched = match filters {
        Some((left, right)) => molframe::analysis::visit_atom_contacts_between(
            &structure,
            &left,
            &right,
            cutoff,
            SpatialBackend::Auto,
            context.execution,
            &mut emit,
        ),
        None => molframe::analysis::visit_atom_contacts(
            &structure,
            cutoff,
            SpatialBackend::Auto,
            context.execution,
            &mut emit,
        ),
    };
    if let Err(error) = searched {
        eprintln!("contact search failed: {error}");
        return Exit::Usage;
    }
    finish_rows(output, write_error.as_ref())
}

fn contact_filters(
    structure: &molframe::Structure,
    between: Option<&[String]>,
    input: &Path,
    context: Context,
) -> Result<Option<(molframe::AtomSelection, molframe::AtomSelection)>, Exit> {
    let Some(between) = between else {
        return Ok(None);
    };
    let [first, second] = between else {
        eprintln!("--between requires exactly two selections");
        return Err(Exit::Usage);
    };
    let first = structure
        .select_text(first, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        })?;
    let second = structure
        .select_text(second, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        })?;
    context.findings(&first.warnings, &input.display().to_string());
    context.findings(&second.warnings, &input.display().to_string());
    Ok(Some((first.selection, second.selection)))
}

pub(crate) fn neighbors(input: &Path, query: &str, cutoff: f32, context: Context) -> Exit {
    use molframe::QueryStructure as _;
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let evaluation = match structure.select_text(query, context.policy, context.execution) {
        Ok(evaluation) => evaluation,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    context.findings(&evaluation.warnings, &input.display().to_string());
    let all = molframe::AtomSelection::All(structure.atom_count());
    let mut output = match RowWriter::new(context, &["atom_a", "atom_b", "distance_angstrom"]) {
        Ok(output) => output,
        Err(error) => return output_error(&error),
    };
    let mut write_error = None;
    let searched = molframe::spatial::for_each_pairs_within_unsorted(
        &molframe::spatial::PairQuery {
            positions: structure.positions(),
            left: &evaluation.selection,
            right: &all,
            cutoff,
            options: molframe::SpatialSearchOptions::with_backend(SpatialBackend::Auto),
            periodic: None,
            context: context.execution,
        },
        |pair| {
            if write_error.is_some() {
                return;
            }
            let first = pair.first.to_string();
            let second = pair.second.to_string();
            let distance = pair.distance_squared.sqrt().to_string();
            write_error = output
                .row([first.as_str(), second.as_str(), distance.as_str()])
                .err();
        },
    );
    if let Err(error) = searched {
        eprintln!("neighbor search failed: {error}");
        return Exit::Usage;
    }
    finish_rows(output, write_error.as_ref())
}

#[derive(Clone, Copy)]
pub(crate) struct SseOptions<'a> {
    pub(crate) ccd: &'a Path,
    pub(crate) ccd_version: &'a str,
    pub(crate) electrostatic_prefactor: f64,
    pub(crate) hydrogen_bond_energy: f64,
    pub(crate) amide_hydrogen_distance: f32,
    pub(crate) minimum_sequence_separation: usize,
    pub(crate) helix_offset: usize,
    pub(crate) turn_offsets: &'a [usize],
}

pub(crate) fn sse(input: &Path, options: SseOptions<'_>, context: Context) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let structure =
        match crate::chemistry::annotate(&structure, options.ccd, options.ccd_version, context) {
            Ok(structure) => structure,
            Err(exit) => return exit,
        };
    let [turn_minimum, turn_maximum] = options.turn_offsets else {
        eprintln!("--turn-offsets requires exactly two values");
        return Exit::Usage;
    };
    let dssp = molframe::analysis::DsspOptions {
        electrostatic_prefactor: options.electrostatic_prefactor,
        hydrogen_bond_energy: options.hydrogen_bond_energy,
        amide_hydrogen_distance: options.amide_hydrogen_distance,
        minimum_sequence_separation: options.minimum_sequence_separation,
        helix_offset: options.helix_offset,
        turn_offsets: *turn_minimum..=*turn_maximum,
    };
    let records = match molframe::analysis::secondary_structure(&structure, &dssp) {
        Ok(records) => records,
        Err(error) => {
            eprintln!("secondary-structure assignment failed: {error}");
            return Exit::Consistency;
        }
    };
    emit_rows(
        context,
        &["residue", "secondary_structure"],
        records,
        |record| {
            vec![
                record.residue.get().to_string(),
                sse_name(record.kind).to_owned(),
            ]
        },
    )
}

pub(crate) fn interfaces(input: &Path, between: &[String], cutoff: f32, context: Context) -> Exit {
    let [first, second] = between else {
        eprintln!("--between requires exactly two chain identifiers");
        return Exit::Usage;
    };
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let residues = match molframe::analysis::chain_interface(
        &structure,
        first,
        second,
        cutoff,
        SpatialBackend::Auto,
        context.execution,
    ) {
        Ok(residues) => residues,
        Err(error) => {
            eprintln!("interface search failed: {error}");
            return Exit::Usage;
        }
    };
    emit_rows(context, &["residue"], residues, |residue| {
        vec![residue.get().to_string()]
    })
}

pub(crate) fn sasa(
    input: &Path,
    probe: f32,
    points: u16,
    radius_set: RadiusSet,
    context: Context,
) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let mut positions = Vec::new();
    let mut radii = Vec::new();
    let mut indices = Vec::new();
    for raw in 0..structure.atom_count() {
        let Some(atom) = structure.atom(AtomIndex::new(raw)) else {
            continue;
        };
        let (Some(position), Some(element)) = (atom.position(), atom.element()) else {
            continue;
        };
        let Some(radius) = molframe::vdw_radius(element, radius_set) else {
            eprintln!("radius set has no value for atom {raw}");
            return Exit::Indeterminate;
        };
        indices.push(raw);
        positions.push(position);
        radii.push(radius);
    }
    let areas = match molframe::surface::shrake_rupley(
        &positions,
        &radii,
        probe,
        points,
        context.execution,
    ) {
        Ok(areas) => areas,
        Err(error) => {
            eprintln!("SASA calculation failed: {error}");
            return Exit::Usage;
        }
    };
    emit_rows(
        context,
        &["atom", "area_angstrom_squared"],
        indices.into_iter().zip(areas),
        |(atom, area)| vec![atom.to_string(), area.to_string()],
    )
}

fn emit_rows<T>(
    context: Context,
    header: &[&str],
    rows: impl IntoIterator<Item = T>,
    render: impl Fn(T) -> Vec<String>,
) -> Exit {
    let mut output = match RowWriter::new(context, header) {
        Ok(output) => output,
        Err(error) => return output_error(&error),
    };
    for item in rows {
        let row = render(item);
        if let Err(error) = output.row(row.iter().map(String::as_str)) {
            return output_error(&error);
        }
    }
    finish_rows(output, None)
}

fn finish_rows(output: RowWriter, error: Option<&std::io::Error>) -> Exit {
    if let Some(error) = error {
        return output_error(error);
    }
    match output.finish() {
        Ok(()) => Exit::Success,
        Err(error) => output_error(&error),
    }
}

fn output_error(error: &std::io::Error) -> Exit {
    eprintln!("could not write result: {error}");
    Exit::Consistency
}

pub(super) fn sse_name(kind: molframe::analysis::SseKind) -> &'static str {
    match kind {
        molframe::analysis::SseKind::AlphaHelix => "alpha-helix",
        molframe::analysis::SseKind::Strand => "strand",
        molframe::analysis::SseKind::Turn => "turn",
        molframe::analysis::SseKind::Coil => "coil",
    }
}
