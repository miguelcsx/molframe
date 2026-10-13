//! What a PDB write cannot carry, reported instead of dropped silently.

use super::PdbOptions;
use super::pdb::{PdbIdentifierNamespace, chain_name};
use molframe_core::bond::BondProvenance;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;

/// One aggregated finding per kind of loss, with a count and the first example.
///
/// Returns nothing when the structure fits the format exactly.
#[must_use]
pub fn projection_loss(structure: &Structure, options: &PdbOptions) -> Vec<Diagnostic> {
    let mut findings = Vec::new();
    let data = structure.data();
    if let Some(id) = data.entry.id.as_deref().filter(|id| id.chars().count() > 4) {
        findings.push(loss(format!(
            "entry identifier `{id}` is longer than four characters and is truncated"
        )));
    }
    let table = &data.topology.entities;
    let described = table
        .iter()
        .filter(|&entity| {
            table.description(entity).is_some() || !table.canonical_sequence(entity).is_empty()
        })
        .count();
    if described > 0 {
        findings.push(loss(format!(
            "{described} entity description(s) and canonical sequence(s) have no PDB \
             counterpart and are not written"
        )));
    }
    renamed_identifiers(structure, options, &mut findings);
    rounded_values(structure, &mut findings);
    unrecorded_bonds(structure, &mut findings);
    findings
}

fn loss(message: String) -> Diagnostic {
    Diagnostic::new(Code::W4101).with_message(message)
}

/// The identifier namespace the options do not write is dropped.
fn renamed_identifiers(
    structure: &Structure,
    options: &PdbOptions,
    findings: &mut Vec<Diagnostic>,
) {
    let namespace = options.identifier_namespace();
    let other = match namespace {
        PdbIdentifierNamespace::Label => PdbIdentifierNamespace::Auth,
        PdbIdentifierNamespace::Auth => PdbIdentifierNamespace::Label,
    };
    let mut differing = 0usize;
    let mut example = None;
    for chain in structure.data().chains() {
        let written = chain_name(chain, namespace);
        let dropped = chain_name(chain, other);
        if dropped.is_some() && dropped != written {
            differing += 1;
            example.get_or_insert((written.unwrap_or(""), dropped.unwrap_or("")));
        }
    }
    if let Some((kept, dropped)) = example {
        findings.push(loss(format!(
            "{differing} chain(s) carry a second {} identifier that is not written \
             (for example `{dropped}` beside written `{kept}`)",
            other.as_str()
        )));
    }
}

#[allow(clippy::cast_possible_truncation)]
fn rounded(value: f32, scale: f64) -> bool {
    let kept = ((f64::from(value) * scale).round() / scale) as f32;
    kept.to_bits() != value.to_bits()
}

fn rounded_values(structure: &Structure, findings: &mut Vec<Diagnostic>) {
    let (mut position, mut occupancy, mut b_factor) = (0usize, 0usize, 0usize);
    for atom in structure.data().atoms() {
        if atom
            .position()
            .is_some_and(|p| p.iter().any(|&v| rounded(v, 1e3)))
        {
            position += 1;
        }
        if atom.occupancy().is_some_and(|v| rounded(v, 1e2)) {
            occupancy += 1;
        }
        if atom.b_factor().is_some_and(|v| rounded(v, 1e2)) {
            b_factor += 1;
        }
    }
    for (count, name, places) in [
        (position, "coordinate", 3),
        (occupancy, "occupancy", 2),
        (b_factor, "B-factor", 2),
    ] {
        if count > 0 {
            findings.push(loss(format!(
                "{count} atom(s) carry a {name} with more than {places} decimals, \
                 rounded in the PDB columns"
            )));
        }
    }
}

fn unrecorded_bonds(structure: &Structure, findings: &mut Vec<Diagnostic>) {
    let count = structure
        .data()
        .bonds
        .iter()
        .filter(|bond| {
            bond.provenance != BondProvenance::File
                && bond.kind.get().is_none_or(|symbol| {
                    structure.resolve(symbol).is_none_or(|kind| {
                        !kind.starts_with("covale")
                            && !matches!(kind, "disulf" | "metalc" | "modres")
                    })
                })
        })
        .count();
    if count > 0 {
        findings.push(loss(format!(
            "{count} bond(s) have no SSBOND, LINK or CONECT record and are not written"
        )));
    }
}

#[cfg(test)]
#[path = "loss_tests.rs"]
mod tests;
