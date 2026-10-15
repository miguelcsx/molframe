//! GROMACS include-topology reader with lossless force-field parameters.
//!
//! Records are parsed once in declaration order. Atom identifiers are validated
//! through an ordered set, giving deterministic duplicate and connectivity
//! errors without requiring identifiers to be contiguous.

use std::collections::{BTreeMap, BTreeSet};

/// Molecule declaration from `[ moleculetype ]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GromacsMoleculeType {
    /// Molecule name.
    pub name: Box<str>,
    /// Number of bonded neighbours excluded from non-bonded interactions.
    pub exclusions: u32,
}

/// One `[ atoms ]` record.
#[derive(Clone, Debug, PartialEq)]
pub struct GromacsItpAtom {
    /// File-level atom number.
    pub number: u32,
    /// Force-field atom type.
    pub atom_type: Box<str>,
    /// Residue number.
    pub residue_number: i64,
    /// Residue name.
    pub residue_name: Box<str>,
    /// Atom name.
    pub atom_name: Box<str>,
    /// Charge-group number.
    pub charge_group: u32,
    /// Partial charge in proton-charge units.
    pub charge: f64,
    /// Optional mass in daltons.
    pub mass: Option<f64>,
    /// Optional free-energy state-B fields, retained verbatim.
    pub state_b: Vec<Box<str>>,
}

/// A typed bonded record with parameters retained as source tokens.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GromacsInteraction<const N: usize> {
    /// Atom numbers in declared order.
    pub atoms: [u32; N],
    /// GROMACS function type.
    pub function: u32,
    /// Remaining force-field parameter tokens.
    pub parameters: Vec<Box<str>>,
}

/// Parsed contents of one GROMACS `.itp` file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GromacsItp {
    /// Preprocessor directives in source order.
    pub directives: Vec<Box<str>>,
    /// Molecule declaration, when present.
    pub molecule_type: Option<GromacsMoleculeType>,
    /// Atom records sorted by atom number.
    pub atoms: Vec<GromacsItpAtom>,
    /// Bond records.
    pub bonds: Vec<GromacsInteraction<2>>,
    /// Pair records.
    pub pairs: Vec<GromacsInteraction<2>>,
    /// Angle records.
    pub angles: Vec<GromacsInteraction<3>>,
    /// Proper and improper dihedral records.
    pub dihedrals: Vec<GromacsInteraction<4>>,
    /// Constraint records.
    pub constraints: Vec<GromacsInteraction<2>>,
    /// Exclusion atom-number lists.
    pub exclusions: Vec<Vec<u32>>,
    /// Unsupported sections retained as comment-free source records.
    pub other_sections: BTreeMap<Box<str>, Vec<Box<str>>>,
}

/// Malformed or inconsistent GROMACS include topology.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum GromacsItpError {
    /// Data appeared before any section header.
    #[error("GROMACS topology data appears outside a section")]
    MissingSection,
    /// A section header was not closed.
    #[error("invalid GROMACS topology section header")]
    InvalidSection,
    /// A record was missing required fields.
    #[error("invalid GROMACS topology record")]
    InvalidRecord,
    /// A numeric field could not be parsed.
    #[error("invalid GROMACS topology numeric field")]
    InvalidNumber,
    /// An atom number was declared more than once.
    #[error("duplicate GROMACS topology atom number")]
    DuplicateAtom,
    /// Bonded or exclusion data references an absent atom number.
    #[error("GROMACS topology interaction references an absent atom")]
    UnknownAtom,
}

/// Parses one already-preprocessed or directive-containing GROMACS `.itp` file.
///
/// Preprocessor directives are retained but not executed. Function-specific
/// parameter tails remain source tokens so the parser does not reinterpret
/// force-field units or discard topology annotations.
///
/// # Errors
///
/// Returns an error for malformed sections or records, duplicate atom numbers,
/// and connectivity to atoms absent from `[ atoms ]`.
pub fn parse_gromacs_itp(text: &str) -> Result<GromacsItp, GromacsItpError> {
    let mut topology = GromacsItp::default();
    let mut sections = BTreeMap::<Box<str>, Vec<Box<str>>>::new();
    let mut current: Option<Box<str>> = None;
    for raw in text.lines() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            topology.directives.push(line.into());
            continue;
        }
        if line.starts_with('[') {
            let Some(end) = line.find(']') else {
                return Err(GromacsItpError::InvalidSection);
            };
            if !line[end + 1..].trim().is_empty() {
                return Err(GromacsItpError::InvalidSection);
            }
            let name: Box<str> = line[1..end].trim().to_ascii_lowercase().into();
            sections.entry(name.clone()).or_default();
            current = Some(name);
            continue;
        }
        let Some(name) = &current else {
            return Err(GromacsItpError::MissingSection);
        };
        let Some(records) = sections.get_mut(name) else {
            return Err(GromacsItpError::MissingSection);
        };
        records.push(line.into());
    }
    populate(&mut topology, sections)?;
    Ok(topology)
}

fn populate(
    topology: &mut GromacsItp,
    mut sections: BTreeMap<Box<str>, Vec<Box<str>>>,
) -> Result<(), GromacsItpError> {
    topology.molecule_type = parse_molecule(sections.remove("moleculetype"))?;
    topology.atoms = parse_atoms(sections.remove("atoms"))?;
    topology.atoms.sort_unstable_by_key(|atom| atom.number);
    let atom_numbers: BTreeSet<_> = topology.atoms.iter().map(|atom| atom.number).collect();
    if atom_numbers.len() != topology.atoms.len() {
        return Err(GromacsItpError::DuplicateAtom);
    }
    topology.bonds = parse_interactions(sections.remove("bonds"), &atom_numbers)?;
    topology.pairs = parse_interactions(sections.remove("pairs"), &atom_numbers)?;
    topology.angles = parse_interactions(sections.remove("angles"), &atom_numbers)?;
    topology.dihedrals = parse_interactions(sections.remove("dihedrals"), &atom_numbers)?;
    topology.constraints = parse_interactions(sections.remove("constraints"), &atom_numbers)?;
    topology.exclusions = parse_exclusions(sections.remove("exclusions"), &atom_numbers)?;
    topology.other_sections = sections;
    Ok(())
}

fn strip_comment(line: &str) -> &str {
    match line.split_once(';') {
        Some((content, _)) => content,
        None => line,
    }
}

fn parse_molecule(
    records: Option<Vec<Box<str>>>,
) -> Result<Option<GromacsMoleculeType>, GromacsItpError> {
    let Some(records) = records else {
        return Ok(None);
    };
    if records.len() != 1 {
        return Err(GromacsItpError::InvalidRecord);
    }
    let fields: Vec<_> = records[0].split_whitespace().collect();
    if fields.len() != 2 {
        return Err(GromacsItpError::InvalidRecord);
    }
    Ok(Some(GromacsMoleculeType {
        name: fields[0].into(),
        exclusions: unsigned(fields[1])?,
    }))
}

fn parse_atoms(records: Option<Vec<Box<str>>>) -> Result<Vec<GromacsItpAtom>, GromacsItpError> {
    let Some(records) = records else {
        return Ok(Vec::new());
    };
    records
        .iter()
        .map(|record| {
            let fields: Vec<_> = record.split_whitespace().collect();
            if fields.len() < 7 {
                return Err(GromacsItpError::InvalidRecord);
            }
            Ok(GromacsItpAtom {
                number: unsigned(fields[0])?,
                atom_type: fields[1].into(),
                residue_number: integer(fields[2])?,
                residue_name: fields[3].into(),
                atom_name: fields[4].into(),
                charge_group: unsigned(fields[5])?,
                charge: number(fields[6])?,
                mass: fields.get(7).map(|value| number(value)).transpose()?,
                state_b: fields.iter().skip(8).map(|value| (*value).into()).collect(),
            })
        })
        .collect()
}

fn parse_interactions<const N: usize>(
    records: Option<Vec<Box<str>>>,
    atom_numbers: &BTreeSet<u32>,
) -> Result<Vec<GromacsInteraction<N>>, GromacsItpError> {
    let Some(records) = records else {
        return Ok(Vec::new());
    };
    records
        .iter()
        .map(|record| {
            let fields: Vec<_> = record.split_whitespace().collect();
            if fields.len() < N + 1 {
                return Err(GromacsItpError::InvalidRecord);
            }
            let mut atoms = [0; N];
            for (index, atom) in atoms.iter_mut().enumerate() {
                *atom = unsigned(fields[index])?;
            }
            if atoms.iter().any(|atom| !atom_numbers.contains(atom)) {
                return Err(GromacsItpError::UnknownAtom);
            }
            Ok(GromacsInteraction {
                atoms,
                function: unsigned(fields[N])?,
                parameters: fields
                    .iter()
                    .skip(N + 1)
                    .map(|value| (*value).into())
                    .collect(),
            })
        })
        .collect()
}

fn parse_exclusions(
    records: Option<Vec<Box<str>>>,
    atom_numbers: &BTreeSet<u32>,
) -> Result<Vec<Vec<u32>>, GromacsItpError> {
    let Some(records) = records else {
        return Ok(Vec::new());
    };
    records
        .iter()
        .map(|record| {
            let values: Vec<_> = record
                .split_whitespace()
                .map(unsigned)
                .collect::<Result<_, _>>()?;
            if values.len() < 2 || values.iter().any(|atom| !atom_numbers.contains(atom)) {
                return Err(GromacsItpError::UnknownAtom);
            }
            Ok(values)
        })
        .collect()
}

fn number(value: &str) -> Result<f64, GromacsItpError> {
    value.parse().map_err(|_| GromacsItpError::InvalidNumber)
}

fn integer(value: &str) -> Result<i64, GromacsItpError> {
    value.parse().map_err(|_| GromacsItpError::InvalidNumber)
}

fn unsigned(value: &str) -> Result<u32, GromacsItpError> {
    value.parse().map_err(|_| GromacsItpError::InvalidNumber)
}

#[cfg(test)]
#[path = "gromacs_itp_tests.rs"]
mod tests;
