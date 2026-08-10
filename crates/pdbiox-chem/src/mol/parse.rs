//! V2000/V3000 parsing and strict multi-record SDF parsing.

use super::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolError, MolRecord, MolVersion, Molecule,
    SdfProperty,
};
use pdbiox_core::element::Element;
use std::collections::BTreeMap;

/// Parses one complete MOL block.
///
/// # Errors
///
/// Returns [`MolError`] when required syntax, numeric fields, or endpoints are invalid.
pub fn parse_mol_record(text: &str) -> Result<MolRecord, MolError> {
    let mut lines = text.lines();
    let name: Box<str> = lines.next().ok_or(MolError::Malformed)?.into();
    let program: Box<str> = lines.next().ok_or(MolError::Malformed)?.into();
    let comment: Box<str> = lines.next().ok_or(MolError::Malformed)?.into();
    let counts = lines.next().ok_or(MolError::Malformed)?;
    if counts.contains("V3000") {
        parse_v3000(name, program, comment, lines)
    } else if counts.contains("V2000") {
        parse_v2000(name, program, comment, counts, lines)
    } else {
        Err(MolError::Malformed)
    }
}

/// Parses every SDF record and refuses a malformed member.
///
/// # Errors
///
/// Returns [`MolError`] for the first malformed member.
pub fn parse_sdf_records(text: &str) -> Result<Vec<MolRecord>, MolError> {
    text.split("$$$$")
        .map(|block| block.trim_start_matches(['\n', '\r']))
        .filter(|block| !block.trim().is_empty())
        .map(parse_mol_record)
        .collect()
}

fn parse_v2000<'a>(
    name: Box<str>,
    program: Box<str>,
    comment: Box<str>,
    counts: &str,
    mut lines: impl Iterator<Item = &'a str>,
) -> Result<MolRecord, MolError> {
    let mut count_fields = counts.split_whitespace();
    let atoms = number::<usize>(count_fields.next())?;
    let bonds = number::<usize>(count_fields.next())?;
    let mut molecule = Molecule {
        atoms: Vec::with_capacity(atoms),
        bonds: Vec::with_capacity(bonds),
    };
    let mut atom_metadata = Vec::with_capacity(atoms);
    for _ in 0..atoms {
        let line = lines.next().ok_or(MolError::Malformed)?;
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(MolError::Malformed);
        }
        molecule.atoms.push(MolAtom {
            position: [
                number(Some(fields[0]))?,
                number(Some(fields[1]))?,
                number(Some(fields[2]))?,
            ],
            element: element(fields[3]),
        });
        atom_metadata.push(MolAtomMetadata {
            formal_charge: fields.get(5).and_then(charge_code),
            stereo_parity: fields.get(6).and_then(nonzero_u8),
            ..MolAtomMetadata::default()
        });
    }
    let mut bond_metadata = Vec::with_capacity(bonds);
    for _ in 0..bonds {
        let fields: Vec<_> = lines
            .next()
            .ok_or(MolError::Malformed)?
            .split_whitespace()
            .collect();
        if fields.len() < 3 {
            return Err(MolError::Malformed);
        }
        let first = one_based(fields[0])?;
        let second = one_based(fields[1])?;
        molecule.bonds.push(MolBond {
            first,
            second,
            order: number(Some(fields[2]))?,
        });
        bond_metadata.push(MolBondMetadata {
            stereo: fields.get(3).and_then(nonzero_u8),
        });
    }
    let mut properties = Vec::new();
    let mut found_end = false;
    while let Some(line) = lines.next() {
        if line == "M  END" {
            found_end = true;
        } else if line.starts_with("M  CHG") && !found_end {
            parse_atom_pairs(line, &mut atom_metadata, |metadata, value| {
                metadata.formal_charge = i8::try_from(value).ok();
            })?;
        } else if line.starts_with("M  ISO") && !found_end {
            parse_atom_pairs(line, &mut atom_metadata, |metadata, value| {
                metadata.isotope = u16::try_from(value).ok();
            })?;
        } else if found_end && let Some(property) = property_header(line) {
            properties.push(parse_property(property, &mut lines));
        }
    }
    if !found_end {
        return Err(MolError::Malformed);
    }
    validate_endpoints(&molecule)?;
    Ok(MolRecord {
        name,
        program,
        comment,
        version: MolVersion::V2000,
        molecule,
        atom_metadata,
        bond_metadata,
        properties,
    })
}

fn parse_v3000<'a>(
    name: Box<str>,
    program: Box<str>,
    comment: Box<str>,
    lines: impl Iterator<Item = &'a str>,
) -> Result<MolRecord, MolError> {
    let mut section = "";
    let mut atoms = Vec::new();
    let mut atom_metadata = Vec::new();
    let mut bonds = Vec::new();
    let mut bond_metadata = Vec::new();
    let mut ids = BTreeMap::new();
    let mut properties = Vec::new();
    let mut counts = None;
    let mut found_ctab_end = false;
    let mut found_m_end = false;
    let mut lines = lines.peekable();
    while let Some(line) = lines.next() {
        let payload = line.strip_prefix("M  V30 ");
        if let Some(payload) = payload {
            if let Some(next) = payload.strip_prefix("BEGIN ") {
                section = next;
                continue;
            }
            if let Some(ended) = payload.strip_prefix("END ") {
                if ended == "CTAB" {
                    found_ctab_end = true;
                }
                section = "";
                continue;
            }
            if let Some(count_line) = payload.strip_prefix("COUNTS ") {
                let fields: Vec<_> = count_line.split_whitespace().collect();
                counts = Some((
                    number::<usize>(fields.first().copied())?,
                    number::<usize>(fields.get(1).copied())?,
                ));
                continue;
            }
            let fields: Vec<_> = payload.split_whitespace().collect();
            if section == "ATOM" {
                if fields.len() < 6 {
                    return Err(MolError::Malformed);
                }
                let id = number::<usize>(fields.first().copied())?;
                if id == 0 || ids.insert(id, atoms.len()).is_some() {
                    return Err(MolError::Malformed);
                }
                atoms.push(MolAtom {
                    element: element(fields[1]),
                    position: [
                        number(Some(fields[2]))?,
                        number(Some(fields[3]))?,
                        number(Some(fields[4]))?,
                    ],
                });
                atom_metadata.push(MolAtomMetadata {
                    formal_charge: tagged(fields.iter().copied(), "CHG=")
                        .and_then(|v| v.parse().ok()),
                    isotope: tagged(fields.iter().copied(), "MASS=").and_then(|v| v.parse().ok()),
                    stereo_parity: tagged(fields.iter().copied(), "CFG=")
                        .and_then(|v| v.parse().ok()),
                });
            } else if section == "BOND" {
                if fields.len() < 4 {
                    return Err(MolError::Malformed);
                }
                let first_id = number::<usize>(Some(fields[2]))?;
                let second_id = number::<usize>(Some(fields[3]))?;
                bonds.push(MolBond {
                    first: *ids.get(&first_id).ok_or(MolError::BondEndpoint)?,
                    second: *ids.get(&second_id).ok_or(MolError::BondEndpoint)?,
                    order: number(Some(fields[1]))?,
                });
                bond_metadata.push(MolBondMetadata {
                    stereo: tagged(fields.iter().copied(), "CFG=").and_then(|v| v.parse().ok()),
                });
            }
        } else if line == "M  END" {
            found_m_end = true;
        } else if found_m_end && let Some(property) = property_header(line) {
            properties.push(parse_property(property, &mut lines));
        }
    }
    if !found_ctab_end || !found_m_end || counts != Some((atoms.len(), bonds.len())) {
        return Err(MolError::Malformed);
    }
    let molecule = Molecule { atoms, bonds };
    validate_endpoints(&molecule)?;
    Ok(MolRecord {
        name,
        program,
        comment,
        version: MolVersion::V3000,
        molecule,
        atom_metadata,
        bond_metadata,
        properties,
    })
}

fn number<T: std::str::FromStr>(value: Option<&str>) -> Result<T, MolError> {
    value
        .and_then(|text| text.parse().ok())
        .ok_or(MolError::Malformed)
}

fn one_based(value: &str) -> Result<usize, MolError> {
    number::<usize>(Some(value))?
        .checked_sub(1)
        .ok_or(MolError::BondEndpoint)
}

fn element(symbol: &str) -> Element {
    match Element::from_symbol(symbol) {
        Some(value) => value,
        None => Element::UNKNOWN,
    }
}

fn nonzero_u8(value: &&str) -> Option<u8> {
    value.parse().ok().filter(|value| *value != 0)
}

fn charge_code(value: &&str) -> Option<i8> {
    match *value {
        "1" => Some(3),
        "2" => Some(2),
        "3" => Some(1),
        "5" => Some(-1),
        "6" => Some(-2),
        "7" => Some(-3),
        _ => None,
    }
}

fn parse_atom_pairs(
    line: &str,
    metadata: &mut [MolAtomMetadata],
    mut assign: impl FnMut(&mut MolAtomMetadata, i32),
) -> Result<(), MolError> {
    let fields: Vec<_> = line.split_whitespace().collect();
    let count = number::<usize>(fields.get(2).copied())?;
    if fields.len() < 3 + count * 2 {
        return Err(MolError::Malformed);
    }
    for pair in fields[3..3 + count * 2].chunks_exact(2) {
        let atom = one_based(pair[0])?;
        let target = metadata.get_mut(atom).ok_or(MolError::BondEndpoint)?;
        assign(target, number(Some(pair[1]))?);
    }
    Ok(())
}

fn property_header(line: &str) -> Option<Box<str>> {
    let start = line.find('<')? + 1;
    let end = line[start..].find('>')? + start;
    Some(line[start..end].into())
}

fn parse_property<'a>(name: Box<str>, lines: &mut impl Iterator<Item = &'a str>) -> SdfProperty {
    let mut value = String::new();
    for line in lines {
        if line.is_empty() {
            break;
        }
        if !value.is_empty() {
            value.push('\n');
        }
        value.push_str(line);
    }
    SdfProperty {
        name,
        value: value.into(),
    }
}

fn tagged<'a>(mut fields: impl Iterator<Item = &'a str>, prefix: &str) -> Option<&'a str> {
    fields.find_map(|field| field.strip_prefix(prefix))
}

fn validate_endpoints(molecule: &Molecule) -> Result<(), MolError> {
    if molecule
        .bonds
        .iter()
        .any(|bond| bond.first >= molecule.atoms.len() || bond.second >= molecule.atoms.len())
    {
        Err(MolError::BondEndpoint)
    } else {
        Ok(())
    }
}
