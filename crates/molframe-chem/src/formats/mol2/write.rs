//! Deterministic Tripos MOL2 writing without synthetic chemistry defaults.

use super::{Mol2Error, Mol2Record};
use std::collections::BTreeSet;
use std::fmt::Write as _;

/// Writes one complete MOL2 molecule record.
///
/// # Errors
///
/// Returns [`Mol2Error`] when fields disagree or cannot be represented without loss.
pub fn write_mol2(record: &Mol2Record) -> Result<String, Mol2Error> {
    validate(record)?;
    let mut output = String::new();
    output.push_str("@<TRIPOS>MOLECULE\n");
    writeln!(output, "{}", record.name).map_err(|_| Mol2Error::Unrepresentable)?;
    write!(
        output,
        "{} {}",
        record.molecule.atoms.len(),
        record.molecule.bonds.len()
    )
    .map_err(|_| Mol2Error::Unrepresentable)?;
    for count in &record.additional_counts {
        write!(output, " {count}").map_err(|_| Mol2Error::Unrepresentable)?;
    }
    output.push('\n');
    writeln!(output, "{}", record.molecule_type).map_err(|_| Mol2Error::Unrepresentable)?;
    writeln!(output, "{}", record.charge_type).map_err(|_| Mol2Error::Unrepresentable)?;
    if let Some(status) = &record.status_bits {
        writeln!(output, "{status}").map_err(|_| Mol2Error::Unrepresentable)?;
    }
    if let Some(comment) = &record.comment {
        if record.status_bits.is_none() {
            output.push('\n');
        }
        writeln!(output, "{comment}").map_err(|_| Mol2Error::Unrepresentable)?;
    }
    output.push_str("@<TRIPOS>ATOM\n");
    for (atom, metadata) in record.molecule.atoms.iter().zip(&record.atom_metadata) {
        write!(
            output,
            "{} {} {:.7} {:.7} {:.7} {}",
            metadata.id,
            metadata.name,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            metadata.atom_type
        )
        .map_err(|_| Mol2Error::Unrepresentable)?;
        write_optional_atom_fields(&mut output, metadata)?;
        output.push('\n');
    }
    output.push_str("@<TRIPOS>BOND\n");
    for (bond, metadata) in record.molecule.bonds.iter().zip(&record.bond_metadata) {
        let first = record.atom_metadata[bond.first].id;
        let second = record.atom_metadata[bond.second].id;
        write!(
            output,
            "{} {} {} {}",
            metadata.id, first, second, metadata.bond_type
        )
        .map_err(|_| Mol2Error::Unrepresentable)?;
        if let Some(status) = &metadata.status_bits {
            write!(output, " {status}").map_err(|_| Mol2Error::Unrepresentable)?;
        }
        output.push('\n');
    }
    for section in &record.extra_sections {
        writeln!(output, "@<TRIPOS>{}", section.name).map_err(|_| Mol2Error::Unrepresentable)?;
        for line in &section.lines {
            writeln!(output, "{line}").map_err(|_| Mol2Error::Unrepresentable)?;
        }
    }
    Ok(output)
}

fn write_optional_atom_fields(
    output: &mut String,
    metadata: &super::Mol2AtomMetadata,
) -> Result<(), Mol2Error> {
    let optional = metadata.substructure_id.is_some()
        || metadata.substructure_name.is_some()
        || metadata.charge.is_some()
        || metadata.status_bits.is_some();
    if !optional {
        return Ok(());
    }
    let substructure_id = metadata.substructure_id.ok_or(Mol2Error::Unrepresentable)?;
    let substructure_name = metadata
        .substructure_name
        .as_deref()
        .ok_or(Mol2Error::Unrepresentable)?;
    write!(output, " {substructure_id} {substructure_name}")
        .map_err(|_| Mol2Error::Unrepresentable)?;
    if metadata.charge.is_some() || metadata.status_bits.is_some() {
        let charge = metadata.charge.ok_or(Mol2Error::Unrepresentable)?;
        write!(output, " {charge:.7}").map_err(|_| Mol2Error::Unrepresentable)?;
    }
    if let Some(status) = &metadata.status_bits {
        write!(output, " {status}").map_err(|_| Mol2Error::Unrepresentable)?;
    }
    Ok(())
}

fn validate(record: &Mol2Record) -> Result<(), Mol2Error> {
    if record.atom_metadata.len() != record.molecule.atoms.len()
        || record.bond_metadata.len() != record.molecule.bonds.len()
    {
        return Err(Mol2Error::MetadataLength);
    }
    let mut atom_ids = BTreeSet::new();
    if record
        .molecule
        .atoms
        .iter()
        .zip(&record.atom_metadata)
        .any(|(atom, metadata)| {
            metadata.id == 0
                || !atom_ids.insert(metadata.id)
                || atom.element == molframe_core::Element::UNKNOWN
                || atom.position.iter().any(|value| !value.is_finite())
                || invalid_token(&metadata.name)
                || invalid_token(&metadata.atom_type)
                || metadata.charge.is_some_and(|value| !value.is_finite())
                || metadata
                    .substructure_name
                    .as_deref()
                    .is_some_and(invalid_token)
        })
    {
        return Err(Mol2Error::Unrepresentable);
    }
    let mut bond_ids = BTreeSet::new();
    if record
        .molecule
        .bonds
        .iter()
        .zip(&record.bond_metadata)
        .any(|(bond, metadata)| {
            metadata.id == 0
                || !bond_ids.insert(metadata.id)
                || bond.first >= record.molecule.atoms.len()
                || bond.second >= record.molecule.atoms.len()
                || invalid_token(&metadata.bond_type)
                || !bond_type_matches(metadata.bond_type.as_ref(), bond.order)
        })
    {
        return Err(Mol2Error::Unrepresentable);
    }
    if invalid_line(&record.name)
        || invalid_token(&record.molecule_type)
        || invalid_token(&record.charge_type)
        || record.status_bits.as_deref().is_some_and(invalid_line)
        || record.comment.as_deref().is_some_and(invalid_line)
        || record.extra_sections.iter().any(|section| {
            invalid_token(&section.name)
                || matches!(section.name.as_ref(), "MOLECULE" | "ATOM" | "BOND")
                || section.lines.iter().any(|line| invalid_line(line))
        })
    {
        return Err(Mol2Error::Unrepresentable);
    }
    Ok(())
}

fn bond_type_matches(kind: &str, order: u8) -> bool {
    match kind {
        "ar" => order == 4,
        "am" => order == 5,
        "du" | "un" | "nc" => order == 0,
        numeric => numeric.parse::<u8>().ok() == Some(order),
    }
}

fn invalid_token(value: &str) -> bool {
    value.is_empty() || value.chars().any(char::is_whitespace)
}

fn invalid_line(value: &str) -> bool {
    value.contains(['\n', '\r'])
}
