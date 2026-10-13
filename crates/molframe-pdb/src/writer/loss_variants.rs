//! What a PQR or PDBQT write cannot carry, reported instead of dropped silently.

use super::PdbOptions;
use super::loss::{loss, renamed_identifiers, rounded, rounded_f64};
use super::pdb::{atom_name, atom_name_field, element_column};
use super::variants::{Variant, present_real, present_symbol};
use molframe_core::SecondaryStructure;
use molframe_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AtomAnnotation, PARTIAL_CHARGE_ANNOTATION,
};
use molframe_core::bond::BondProvenance;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::element::Element;
use molframe_core::structure::Structure;

/// One aggregated finding per kind of loss a PQR write incurs.
#[must_use]
pub fn pqr_loss(structure: &Structure, options: &PdbOptions) -> Vec<Diagnostic> {
    variant_loss(structure, options, Variant::Pqr)
}

/// One aggregated finding per kind of loss a rigid PDBQT write incurs.
#[must_use]
pub fn pdbqt_loss(structure: &Structure, options: &PdbOptions) -> Vec<Diagnostic> {
    variant_loss(structure, options, Variant::Pdbqt)
}

fn variant_loss(structure: &Structure, options: &PdbOptions, variant: Variant) -> Vec<Diagnostic> {
    let mut findings = Vec::new();
    let name = if variant == Variant::Pqr {
        "PQR"
    } else {
        "PDBQT"
    };
    let data = structure.data();
    let entry = &data.entry;
    if entry.id.is_some()
        || entry.title.is_some()
        || entry.method.is_some()
        || entry.resolution.is_some()
        || data.cell.is_some()
    {
        findings.push(loss(format!(
            "the entry identifier, title, method, resolution and cell have no {name} record \
             and are not written"
        )));
    }
    if data
        .bonds
        .iter()
        .any(|bond| bond.provenance == BondProvenance::File)
        || structure
            .secondary_structure()
            .iter()
            .any(|&class| class != SecondaryStructure::Unknown)
    {
        findings.push(loss(format!(
            "bonds and secondary structure have no {name} record and are not written"
        )));
    }
    renamed_identifiers(structure, options, &mut findings);
    atom_losses(structure, options, variant, name, &mut findings);
    findings
}

fn atom_losses(
    structure: &Structure,
    options: &PdbOptions,
    variant: Variant,
    name: &str,
    findings: &mut Vec<Diagnostic>,
) {
    let charges = real_column(structure, PARTIAL_CHARGE_ANNOTATION);
    let radii = real_column(structure, ATOM_RADIUS_ANNOTATION);
    let types = match structure.annotations().get(AUTODOCK_TYPE_ANNOTATION) {
        Some(AtomAnnotation::Symbol(column)) => Some(column),
        _ => None,
    };
    let charge_places = if variant == Variant::Pqr { 4 } else { 3 };
    let counts = &mut [0usize; 7];
    for atom in structure.data().atoms() {
        let index = atom.index().get();
        if atom.formal_charge().is_some_and(|charge| charge != 0) {
            counts[0] += 1;
        }
        if atom
            .position()
            .is_some_and(|p| p.iter().any(|&v| rounded(v, 1e3)))
        {
            counts[1] += 1;
        }
        let scale = 10f64.powi(charge_places);
        if charges
            .and_then(|column| present_real(column, index))
            .is_some_and(|value| rounded_f64(value, scale))
        {
            counts[2] += 1;
        }
        if variant == Variant::Pqr {
            if radii
                .and_then(|column| present_real(column, index))
                .is_some_and(|value| rounded_f64(value, 1e4))
            {
                counts[3] += 1;
            }
            if atom
                .occupancy()
                .is_some_and(|v| v.to_bits() != 1.0f32.to_bits())
                || atom.b_factor().is_some_and(|v| v != 0.0)
            {
                counts[4] += 1;
            }
        } else {
            if atom.occupancy().is_some_and(|v| rounded(v, 1e2))
                || atom.b_factor().is_some_and(|v| rounded(v, 1e2))
            {
                counts[3] += 1;
            }
            if types
                .and_then(|column| present_symbol(column, index))
                .and_then(|symbol| structure.resolve(symbol))
                .is_some_and(|symbol| symbol.chars().count() > 2)
            {
                counts[5] += 1;
            }
        }
        if let (Some(element), Some(label)) = (
            atom.element(),
            atom_name(atom, options.identifier_namespace()),
        ) {
            let field = atom_name_field(label, &element_column(&atom));
            if element != Element::UNKNOWN && Element::infer_from_pdb_atom_name(&field) != element {
                counts[6] += 1;
            }
        }
    }
    report(counts, variant, name, charge_places, findings);
}

fn report(
    counts: &[usize; 7],
    variant: Variant,
    name: &str,
    charge_places: i32,
    findings: &mut Vec<Diagnostic>,
) {
    let messages = [
        (
            counts[0],
            "a non-zero formal charge, which has no column".to_owned(),
        ),
        (
            counts[1],
            "a coordinate with more than 3 decimals, rounded".to_owned(),
        ),
        (
            counts[2],
            format!("a partial charge with more than {charge_places} decimals, rounded"),
        ),
        (
            counts[3],
            if variant == Variant::Pqr {
                "a radius with more than 4 decimals, rounded".to_owned()
            } else {
                "an occupancy or B-factor with more than 2 decimals, rounded".to_owned()
            },
        ),
        (
            counts[4],
            "an occupancy or B-factor that PQR has no column for".to_owned(),
        ),
        (
            counts[5],
            "an AutoDock type longer than its two columns, which shifts the row".to_owned(),
        ),
        (
            counts[6],
            "an element its atom name does not imply, so it reads back inferred".to_owned(),
        ),
    ];
    for (count, text) in messages {
        if count > 0 {
            findings.push(loss(format!("{count} atom(s) carry {text} in {name}")));
        }
    }
}

fn real_column<'a>(
    structure: &'a Structure,
    name: &str,
) -> Option<&'a molframe_core::annotation::AnnotationColumn<f64>> {
    match structure.annotations().get(name) {
        Some(AtomAnnotation::Real(column)) => Some(column),
        _ => None,
    }
}

#[cfg(test)]
#[path = "loss_variants_tests.rs"]
mod tests;
