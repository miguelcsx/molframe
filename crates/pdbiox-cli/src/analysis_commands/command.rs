//! Thin rendering over native analysis kernels.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use pdbiox::QueryStructure as _;
use pdbiox::{AtomIndex, RadiusSet, SpatialBackend};
use std::fmt::Write as _;
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
    let contacts = match pdbiox::analysis::atom_contacts(&structure, cutoff, SpatialBackend::Auto) {
        Ok(contacts) => contacts,
        Err(error) => {
            eprintln!("contact search failed: {error}");
            return Exit::Usage;
        }
    };
    let filters = match contact_filters(&structure, between, input, context) {
        Ok(filters) => filters,
        Err(exit) => return exit,
    };
    let rows = contacts
        .iter()
        .filter(|contact| contact_allowed(contact.first, contact.second, filters.as_ref()))
        .map(|contact| {
            vec![
                contact.first.get().to_string(),
                contact.second.get().to_string(),
                contact.distance.to_string(),
            ]
        })
        .collect::<Vec<_>>();
    emit_rows(
        context,
        &["first_atom", "second_atom", "distance_angstrom"],
        &rows,
    );
    Exit::Success
}

fn contact_filters(
    structure: &pdbiox::Structure,
    between: Option<&[String]>,
    input: &Path,
    context: Context,
) -> Result<Option<(pdbiox::AtomSelection, pdbiox::AtomSelection)>, Exit> {
    let Some(between) = between else {
        return Ok(None);
    };
    let [first, second] = between else {
        eprintln!("--between requires exactly two selections");
        return Err(Exit::Usage);
    };
    let first = structure
        .select_text(first, context.policy)
        .map_err(|findings| {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        })?;
    let second = structure
        .select_text(second, context.policy)
        .map_err(|findings| {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        })?;
    context.findings(&first.warnings, &input.display().to_string());
    context.findings(&second.warnings, &input.display().to_string());
    Ok(Some((first.selection, second.selection)))
}

fn contact_allowed(
    first: AtomIndex,
    second: AtomIndex,
    filters: Option<&(pdbiox::AtomSelection, pdbiox::AtomSelection)>,
) -> bool {
    filters.is_none_or(|(left, right)| {
        (left.contains(first.get()) && right.contains(second.get()))
            || (left.contains(second.get()) && right.contains(first.get()))
    })
}

pub(crate) fn neighbors(input: &Path, query: &str, cutoff: f32, context: Context) -> Exit {
    use pdbiox::QueryStructure as _;
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let evaluation = match structure.select_text(query, context.policy) {
        Ok(evaluation) => evaluation,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    context.findings(&evaluation.warnings, &input.display().to_string());
    let all = pdbiox::AtomSelection::All(structure.atom_count());
    let pairs = match pdbiox::pairs_within(
        structure.positions(),
        &evaluation.selection,
        &all,
        cutoff,
        SpatialBackend::Auto,
        None,
    ) {
        Ok(pairs) => pairs,
        Err(error) => {
            eprintln!("neighbor search failed: {error}");
            return Exit::Usage;
        }
    };
    let rows = pairs
        .into_iter()
        .map(|pair| {
            vec![
                pair.first.to_string(),
                pair.second.to_string(),
                pair.distance_squared.sqrt().to_string(),
            ]
        })
        .collect::<Vec<_>>();
    emit_rows(context, &["atom_a", "atom_b", "distance_angstrom"], &rows);
    Exit::Success
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
    let dssp = pdbiox::analysis::DsspOptions {
        electrostatic_prefactor: options.electrostatic_prefactor,
        hydrogen_bond_energy: options.hydrogen_bond_energy,
        amide_hydrogen_distance: options.amide_hydrogen_distance,
        minimum_sequence_separation: options.minimum_sequence_separation,
        helix_offset: options.helix_offset,
        turn_offsets: *turn_minimum..=*turn_maximum,
    };
    let records = match pdbiox::analysis::secondary_structure(&structure, &dssp) {
        Ok(records) => records,
        Err(error) => {
            eprintln!("secondary-structure assignment failed: {error}");
            return Exit::Consistency;
        }
    };
    let rows = records
        .into_iter()
        .map(|record| {
            vec![
                record.residue.get().to_string(),
                sse_name(record.kind).to_owned(),
            ]
        })
        .collect::<Vec<_>>();
    emit_rows(context, &["residue", "secondary_structure"], &rows);
    Exit::Success
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
    let residues = match pdbiox::analysis::chain_interface(
        &structure,
        first,
        second,
        cutoff,
        SpatialBackend::Auto,
    ) {
        Ok(residues) => residues,
        Err(error) => {
            eprintln!("interface search failed: {error}");
            return Exit::Usage;
        }
    };
    let rows = residues
        .into_iter()
        .map(|residue| vec![residue.get().to_string()])
        .collect::<Vec<_>>();
    emit_rows(context, &["residue"], &rows);
    Exit::Success
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
        let Some(radius) = pdbiox::vdw_radius(element, radius_set) else {
            eprintln!("radius set has no value for atom {raw}");
            return Exit::Indeterminate;
        };
        indices.push(raw);
        positions.push(position);
        radii.push(radius);
    }
    let areas = match pdbiox::surface::shrake_rupley(&positions, &radii, probe, points) {
        Ok(areas) => areas,
        Err(error) => {
            eprintln!("SASA calculation failed: {error}");
            return Exit::Usage;
        }
    };
    let rows = indices
        .into_iter()
        .zip(areas)
        .map(|(atom, area)| vec![atom.to_string(), area.to_string()])
        .collect::<Vec<_>>();
    emit_rows(context, &["atom", "area_angstrom_squared"], &rows);
    Exit::Success
}

fn emit_rows(context: Context, header: &[&str], rows: &[Vec<String>]) {
    if context.is_json() {
        let objects = rows
            .iter()
            .map(|row| {
                let mut json = Json::new();
                for (key, value) in header.iter().zip(row) {
                    json.text(key, value);
                }
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, header);
        for row in rows {
            table.row(row.iter().map(String::as_str));
        }
        context.result(&table.finish());
    } else {
        let mut text = String::new();
        let _ = writeln!(text, "{}", header.join("\t"));
        for row in rows {
            let _ = writeln!(text, "{}", row.join("\t"));
        }
        context.result(text.trim_end());
    }
}

pub(super) fn sse_name(kind: pdbiox::analysis::SseKind) -> &'static str {
    match kind {
        pdbiox::analysis::SseKind::AlphaHelix => "alpha-helix",
        pdbiox::analysis::SseKind::Strand => "strand",
        pdbiox::analysis::SseKind::Turn => "turn",
        pdbiox::analysis::SseKind::Coil => "coil",
    }
}
