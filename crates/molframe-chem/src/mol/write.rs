//! Deterministic V2000/V3000 and SDF writing.

use super::{MolError, MolRecord, MolVersion, SdfProperty};
use std::fmt::Write as _;

/// Writes one MOL block in its declared dialect.
///
/// # Errors
///
/// Returns [`MolError`] when the record is inconsistent or its selected dialect
/// cannot represent a value.
pub fn write_mol(record: &MolRecord) -> Result<String, MolError> {
    validate(record)?;
    let mut output = String::new();
    writeln!(output, "{}", record.name).map_err(|_| MolError::Malformed)?;
    writeln!(output, "{}", record.program).map_err(|_| MolError::Malformed)?;
    writeln!(output, "{}", record.comment).map_err(|_| MolError::Malformed)?;
    match record.version {
        MolVersion::V2000 => write_v2000(record, &mut output)?,
        MolVersion::V3000 => write_v3000(record, &mut output)?,
    }
    write_properties(&record.properties, &mut output)?;
    Ok(output)
}

/// Writes ordered SDF records with an explicit terminator per record.
///
/// # Errors
///
/// Returns [`MolError`] when any member cannot be written faithfully.
pub fn write_sdf(records: &[MolRecord]) -> Result<String, MolError> {
    let mut output = String::new();
    for record in records {
        output.push_str(&write_mol(record)?);
        output.push_str("$$$$\n");
    }
    Ok(output)
}

fn validate(record: &MolRecord) -> Result<(), MolError> {
    if record.atom_metadata.len() != record.molecule.atoms.len()
        || record.bond_metadata.len() != record.molecule.bonds.len()
    {
        return Err(MolError::MetadataLength);
    }
    if record.molecule.atoms.iter().any(|atom| {
        atom.element == molframe_core::element::Element::UNKNOWN
            || atom.position.iter().any(|value| !value.is_finite())
    }) {
        return Err(MolError::Unrepresentable(record.version));
    }
    if record.molecule.bonds.iter().any(|bond| {
        bond.first >= record.molecule.atoms.len()
            || bond.second >= record.molecule.atoms.len()
            || !(1..=4).contains(&bond.order)
    }) {
        return Err(MolError::BondEndpoint);
    }
    if [&record.name, &record.program, &record.comment]
        .into_iter()
        .any(|value| value.contains(['\n', '\r']))
        || record.properties.iter().any(invalid_property)
    {
        return Err(MolError::Unrepresentable(record.version));
    }
    Ok(())
}

fn write_v2000(record: &MolRecord, output: &mut String) -> Result<(), MolError> {
    if record.molecule.atoms.len() > 999 || record.molecule.bonds.len() > 999 {
        return Err(MolError::Unrepresentable(MolVersion::V2000));
    }
    if record
        .atom_metadata
        .iter()
        .any(|metadata| metadata.isotope.is_some_and(|value| value > 9_999))
    {
        return Err(MolError::Unrepresentable(MolVersion::V2000));
    }
    writeln!(
        output,
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
        record.molecule.atoms.len(),
        record.molecule.bonds.len()
    )
    .map_err(|_| MolError::Malformed)?;
    for (atom, metadata) in record.molecule.atoms.iter().zip(&record.atom_metadata) {
        let x = coordinate(atom.position[0], MolVersion::V2000)?;
        let y = coordinate(atom.position[1], MolVersion::V2000)?;
        let z = coordinate(atom.position[2], MolVersion::V2000)?;
        let parity = match metadata.stereo_parity {
            Some(value) => value,
            None => 0,
        };
        writeln!(
            output,
            "{x}{y}{z} {:<3} 0  0 {:>2}  0  0  0  0  0  0  0  0  0  0",
            atom.element.symbol(),
            parity
        )
        .map_err(|_| MolError::Malformed)?;
    }
    for (bond, metadata) in record.molecule.bonds.iter().zip(&record.bond_metadata) {
        let stereo = match metadata.stereo {
            Some(value) => value,
            None => 0,
        };
        writeln!(
            output,
            "{:>3}{:>3}{:>3}{:>3}  0  0  0",
            bond.first + 1,
            bond.second + 1,
            bond.order,
            stereo
        )
        .map_err(|_| MolError::Malformed)?;
    }
    write_atom_pairs(record, output, "CHG", |metadata| {
        metadata.formal_charge.map(i32::from)
    })?;
    write_atom_pairs(record, output, "ISO", |metadata| {
        metadata.isotope.map(i32::from)
    })?;
    output.push_str("M  END\n");
    Ok(())
}

fn write_v3000(record: &MolRecord, output: &mut String) -> Result<(), MolError> {
    output.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n");
    writeln!(
        output,
        "M  V30 COUNTS {} {} 0 0 0",
        record.molecule.atoms.len(),
        record.molecule.bonds.len()
    )
    .map_err(|_| MolError::Malformed)?;
    output.push_str("M  V30 BEGIN ATOM\n");
    for (index, (atom, metadata)) in record
        .molecule
        .atoms
        .iter()
        .zip(&record.atom_metadata)
        .enumerate()
    {
        write!(
            output,
            "M  V30 {} {} {:.7} {:.7} {:.7} 0",
            index + 1,
            atom.element.symbol(),
            atom.position[0],
            atom.position[1],
            atom.position[2]
        )
        .map_err(|_| MolError::Malformed)?;
        if let Some(charge) = metadata.formal_charge {
            write!(output, " CHG={charge}").map_err(|_| MolError::Malformed)?;
        }
        if let Some(isotope) = metadata.isotope {
            write!(output, " MASS={isotope}").map_err(|_| MolError::Malformed)?;
        }
        if let Some(parity) = metadata.stereo_parity {
            write!(output, " CFG={parity}").map_err(|_| MolError::Malformed)?;
        }
        output.push('\n');
    }
    output.push_str("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for (index, (bond, metadata)) in record
        .molecule
        .bonds
        .iter()
        .zip(&record.bond_metadata)
        .enumerate()
    {
        write!(
            output,
            "M  V30 {} {} {} {}",
            index + 1,
            bond.order,
            bond.first + 1,
            bond.second + 1
        )
        .map_err(|_| MolError::Malformed)?;
        if let Some(stereo) = metadata.stereo {
            write!(output, " CFG={stereo}").map_err(|_| MolError::Malformed)?;
        }
        output.push('\n');
    }
    output.push_str("M  V30 END BOND\nM  V30 END CTAB\nM  END\n");
    Ok(())
}

fn coordinate(value: f32, version: MolVersion) -> Result<String, MolError> {
    let rendered = format!("{value:>10.4}");
    if rendered.len() > 10 {
        Err(MolError::Unrepresentable(version))
    } else {
        Ok(rendered)
    }
}

fn write_atom_pairs(
    record: &MolRecord,
    output: &mut String,
    kind: &str,
    value: impl Fn(&super::MolAtomMetadata) -> Option<i32>,
) -> Result<(), MolError> {
    let pairs: Vec<_> = record
        .atom_metadata
        .iter()
        .enumerate()
        .filter_map(|(index, metadata)| value(metadata).map(|value| (index + 1, value)))
        .collect();
    for chunk in pairs.chunks(8) {
        write!(output, "M  {kind}{:>3}", chunk.len()).map_err(|_| MolError::Malformed)?;
        for (atom, value) in chunk {
            write!(output, "{atom:>4}{value:>4}").map_err(|_| MolError::Malformed)?;
        }
        output.push('\n');
    }
    Ok(())
}

fn write_properties(properties: &[SdfProperty], output: &mut String) -> Result<(), MolError> {
    for property in properties {
        if property.name.is_empty() || property.name.contains(['<', '>', '\n', '\r']) {
            return Err(MolError::Malformed);
        }
        writeln!(output, ">  <{}>", property.name).map_err(|_| MolError::Malformed)?;
        writeln!(output, "{}\n", property.value).map_err(|_| MolError::Malformed)?;
    }
    Ok(())
}

fn invalid_property(property: &SdfProperty) -> bool {
    property.name.is_empty()
        || property.name.contains(['<', '>', '\n', '\r'])
        || property.value.contains('\r')
        || property.value.starts_with('\n')
        || property.value.ends_with('\n')
        || property.value.contains("\n\n")
        || property.value.lines().any(|line| line == "$$$$")
}
