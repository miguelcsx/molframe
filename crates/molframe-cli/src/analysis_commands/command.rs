//! Thin rendering over native analysis kernels.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, RowWriter, emit_rows, finish_rows, output_error};
use molframe::chemistry::RadiusSet;
use molframe::spatial::SpatialBackend;
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
            structure.engine(),
            &left,
            &right,
            cutoff,
            SpatialBackend::Auto,
            context.execution,
            &mut emit,
        ),
        None => molframe::analysis::visit_atom_contacts(
            structure.engine(),
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
) -> Result<
    Option<(
        molframe_core::selection::AtomSelection,
        molframe_core::selection::AtomSelection,
    )>,
    Exit,
> {
    let Some(between) = between else {
        return Ok(None);
    };
    let [first, second] = between else {
        eprintln!("--between requires exactly two selections");
        return Err(Exit::Usage);
    };
    let first = crate::commands::select_text(structure, first, context.policy, context.execution)
        .map_err(|findings| {
        context.findings(&findings, &input.display().to_string());
        Exit::of(&findings)
    })?;
    let second = crate::commands::select_text(structure, second, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        })?;
    context.findings(&first.warnings, &input.display().to_string());
    context.findings(&second.warnings, &input.display().to_string());
    Ok(Some((first.selection, second.selection)))
}

pub(crate) fn neighbors(input: &Path, query: &str, cutoff: f32, context: Context) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let evaluation =
        match crate::commands::select_text(&structure, query, context.policy, context.execution) {
            Ok(evaluation) => evaluation,
            Err(findings) => {
                context.findings(&findings, &input.display().to_string());
                return Exit::of(&findings);
            }
        };
    context.findings(&evaluation.warnings, &input.display().to_string());
    let all = molframe_core::selection::AtomSelection::All(structure.atom_count());
    let mut output = match RowWriter::new(context, &["atom_a", "atom_b", "distance_angstrom"]) {
        Ok(output) => output,
        Err(error) => return output_error(&error),
    };
    let mut write_error = None;
    let searched = molframe::spatial::for_each_pairs_within_unsorted(
        &molframe::spatial::PairQuery {
            positions: structure.coordinates(),
            left: &evaluation.selection,
            right: &all,
            cutoff,
            options: molframe::spatial::SpatialSearchOptions::with_backend(SpatialBackend::Auto),
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
    pub(crate) role_profile: &'a Path,
    pub(crate) electrostatic_prefactor: f64,
    pub(crate) hydrogen_bond_energy: f64,
    pub(crate) amide_hydrogen_distance: f32,
    pub(crate) minimum_sequence_separation: usize,
    pub(crate) helix_offset: usize,
    pub(crate) three_ten_offset: usize,
    pub(crate) pi_offset: usize,
    pub(crate) turn_offsets: &'a [usize],
    pub(crate) bend_angle_degrees: f32,
}

pub(crate) fn sse(input: &Path, options: SseOptions<'_>, context: Context) -> Exit {
    let profile = match crate::role_profile::read(options.role_profile) {
        Ok(profile) => profile,
        Err(error) => {
            eprintln!("polymer role profile failed: {error}");
            return Exit::Policy;
        }
    };
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let provider = match crate::chemistry::load_ccd(options.ccd, options.ccd_version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let roles = match molframe::chemistry::apply_polymer_role_profile(
        structure.engine(),
        &provider,
        &profile,
    ) {
        Ok(report) => report,
        Err(finding) => {
            context.findings(&[finding], &options.role_profile.display().to_string());
            return Exit::Policy;
        }
    };
    for component in &roles.unresolved_components {
        eprintln!(
            "polymer role profile {}: unresolved CCD component {component}",
            roles.profile_id
        );
    }
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
        three_ten_offset: options.three_ten_offset,
        pi_offset: options.pi_offset,
        turn_offsets: *turn_minimum..=*turn_maximum,
        bend_angle_degrees: options.bend_angle_degrees,
    };
    let records = match molframe::analysis::secondary_structure(&roles.structure, &dssp) {
        Ok(records) => records,
        Err(error) => {
            eprintln!("secondary-structure assignment failed: {error}");
            return Exit::Consistency;
        }
    };
    emit_rows(
        context,
        &["residue", "secondary_structure"],
        records.iter(),
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
        structure.engine(),
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
        let Some(atom) = structure.atom_at(raw as usize) else {
            continue;
        };
        let (Some(position), Some(element)) = (atom.position(), atom.element()) else {
            continue;
        };
        let Some(radius) = molframe::chemistry::vdw_radius(element, radius_set) else {
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

pub(super) const fn sse_name(kind: molframe::SecondaryStructure) -> &'static str {
    use molframe::SecondaryStructure as Ss;
    match kind {
        Ss::Unknown => "unknown",
        Ss::Coil => "coil",
        Ss::AlphaHelix => "alpha-helix",
        Ss::ThreeTenHelix => "3-10-helix",
        Ss::PiHelix => "pi-helix",
        Ss::PolyProline => "polyproline",
        Ss::OtherHelix => "other-helix",
        Ss::BetaBridge => "beta-bridge",
        Ss::Strand => "strand",
        Ss::Turn => "turn",
        Ss::Bend => "bend",
    }
}
