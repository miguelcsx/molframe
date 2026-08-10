//! Strict Tripos MOL2 parsing with source identifiers and extra sections.

use super::{Mol2AtomMetadata, Mol2BondMetadata, Mol2Error, Mol2Record, Mol2Section};
use crate::{MolAtom, MolBond, Molecule};
use pdbiox_core::element::Element;
use std::collections::{BTreeMap, BTreeSet};

const AROMATIC_ORDER: u8 = 4;
const AMIDE_ORDER: u8 = 5;

type ParsedAtoms = (Molecule, Vec<Mol2AtomMetadata>, BTreeMap<usize, usize>);

/// Parses one complete MOL2 molecule record.
///
/// # Errors
///
/// Returns [`Mol2Error`] for malformed syntax, inconsistent counts, or invalid identifiers.
pub fn parse_mol2_record(text: &str) -> Result<Mol2Record, Mol2Error> {
    let sections = sections(text)?;
    let molecule_lines = required_section(&sections, "MOLECULE")?;
    if molecule_lines.len() < 4 {
        return Err(Mol2Error::Malformed);
    }
    let counts: Vec<_> = molecule_lines[1].split_whitespace().collect();
    let atom_count = number(counts.first().copied())?;
    let bond_count = number(counts.get(1).copied())?;
    let additional_counts = counts[2..]
        .iter()
        .map(|value| number(Some(value)))
        .collect::<Result<Vec<_>, _>>()?;
    let name = molecule_lines[0].clone();
    let molecule_type = molecule_lines[2].clone();
    let charge_type = molecule_lines[3].clone();
    let status_bits = molecule_lines
        .get(4)
        .cloned()
        .filter(|value| !value.is_empty());
    let comment = molecule_lines
        .get(5)
        .cloned()
        .filter(|value| !value.is_empty());
    let (molecule, atom_metadata, atom_ids) = parse_atoms(required_section(&sections, "ATOM")?)?;
    let (bonds, bond_metadata) = parse_bonds(required_section(&sections, "BOND")?, &atom_ids)?;
    if molecule.atoms.len() != atom_count || bonds.len() != bond_count {
        return Err(Mol2Error::CountMismatch);
    }
    let extra_sections = sections
        .into_iter()
        .filter(|section| !matches!(section.name.as_ref(), "MOLECULE" | "ATOM" | "BOND"))
        .collect();
    Ok(Mol2Record {
        name,
        molecule_type,
        charge_type,
        additional_counts,
        status_bits,
        comment,
        molecule: Molecule {
            atoms: molecule.atoms,
            bonds,
        },
        atom_metadata,
        bond_metadata,
        extra_sections,
    })
}

fn sections(text: &str) -> Result<Vec<Mol2Section>, Mol2Error> {
    let mut result: Vec<Mol2Section> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.trim().strip_prefix("@<TRIPOS>") {
            if name.is_empty() || result.iter().any(|section| section.name.as_ref() == name) {
                return Err(Mol2Error::Malformed);
            }
            result.push(Mol2Section {
                name: name.into(),
                lines: Vec::new(),
            });
        } else if let Some(section) = result.last_mut() {
            section.lines.push(line.into());
        } else if !line.trim().is_empty() {
            return Err(Mol2Error::Malformed);
        }
    }
    Ok(result)
}

fn required_section<'a>(
    sections: &'a [Mol2Section],
    name: &str,
) -> Result<&'a [Box<str>], Mol2Error> {
    sections
        .iter()
        .find(|section| section.name.as_ref() == name)
        .map(|section| section.lines.as_slice())
        .ok_or(Mol2Error::Malformed)
}

fn parse_atoms(lines: &[Box<str>]) -> Result<ParsedAtoms, Mol2Error> {
    let mut molecule = Molecule::default();
    let mut metadata = Vec::new();
    let mut ids = BTreeMap::new();
    for line in lines.iter().filter(|line| !line.trim().is_empty()) {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(Mol2Error::Malformed);
        }
        let id = positive(fields[0])?;
        if ids.insert(id, molecule.atoms.len()).is_some() {
            return Err(Mol2Error::InvalidIdentifier);
        }
        let position = [
            number(Some(fields[2]))?,
            number(Some(fields[3]))?,
            number(Some(fields[4]))?,
        ];
        molecule.atoms.push(MolAtom {
            element: element(fields[5]),
            position,
        });
        metadata.push(Mol2AtomMetadata {
            id,
            name: fields[1].into(),
            atom_type: fields[5].into(),
            substructure_id: fields.get(6).map(|value| positive(value)).transpose()?,
            substructure_name: fields.get(7).map(|value| (*value).into()),
            charge: fields.get(8).map(|value| number(Some(value))).transpose()?,
            status_bits: fields
                .get(9..)
                .filter(|tail| !tail.is_empty())
                .map(|tail| tail.join(" ").into()),
        });
    }
    Ok((molecule, metadata, ids))
}

fn parse_bonds(
    lines: &[Box<str>],
    atom_ids: &BTreeMap<usize, usize>,
) -> Result<(Vec<MolBond>, Vec<Mol2BondMetadata>), Mol2Error> {
    let mut bonds = Vec::new();
    let mut metadata = Vec::new();
    let mut ids = BTreeSet::new();
    for line in lines.iter().filter(|line| !line.trim().is_empty()) {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(Mol2Error::Malformed);
        }
        let id = positive(fields[0])?;
        if !ids.insert(id) {
            return Err(Mol2Error::InvalidIdentifier);
        }
        let first = atom_ids
            .get(&positive(fields[1])?)
            .copied()
            .ok_or(Mol2Error::InvalidIdentifier)?;
        let second = atom_ids
            .get(&positive(fields[2])?)
            .copied()
            .ok_or(Mol2Error::InvalidIdentifier)?;
        bonds.push(MolBond {
            first,
            second,
            order: bond_order(fields[3])?,
        });
        metadata.push(Mol2BondMetadata {
            id,
            bond_type: fields[3].into(),
            status_bits: fields
                .get(4..)
                .filter(|tail| !tail.is_empty())
                .map(|tail| tail.join(" ").into()),
        });
    }
    Ok((bonds, metadata))
}

fn bond_order(value: &str) -> Result<u8, Mol2Error> {
    match value {
        "ar" => Ok(AROMATIC_ORDER),
        "am" => Ok(AMIDE_ORDER),
        "du" | "un" | "nc" => Ok(0),
        value => value.parse().map_err(|_| Mol2Error::Malformed),
    }
}

fn element(atom_type: &str) -> Element {
    let symbol = atom_type
        .split_once('.')
        .map_or(atom_type, |(symbol, _)| symbol);
    match Element::from_symbol(symbol) {
        Some(element) => element,
        None => Element::UNKNOWN,
    }
}

fn positive(value: &str) -> Result<usize, Mol2Error> {
    number(Some(value)).and_then(|value| {
        if value == 0 {
            Err(Mol2Error::InvalidIdentifier)
        } else {
            Ok(value)
        }
    })
}

fn number<T: std::str::FromStr>(value: Option<&str>) -> Result<T, Mol2Error> {
    value
        .and_then(|text| text.parse().ok())
        .ok_or(Mol2Error::Malformed)
}
