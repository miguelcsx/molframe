//! LAMMPS text data topology reader.
//!
//! Atom identifiers remain the file's identifiers because LAMMPS permits sparse
//! numbering and connectivity refers to those values. Parsing is linear in the
//! number of records; identifier validation uses an ordered set for deterministic
//! duplicate and reference checks.

use std::collections::{BTreeMap, BTreeSet};

/// Atom style declared on the `Atoms` section.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum LammpsAtomStyle {
    /// `id type x y z`.
    Atomic,
    /// `id type charge x y z`.
    Charge,
    /// `id molecule type x y z`.
    Molecular,
    /// `id molecule type charge x y z`.
    Full,
}

/// One atom in a LAMMPS data topology.
#[derive(Clone, Debug, PartialEq)]
pub struct LammpsDataAtom {
    /// File-level atom identifier.
    pub id: i64,
    /// Optional molecule identifier.
    pub molecule: Option<i64>,
    /// Force-field atom type.
    pub atom_type: u32,
    /// Optional partial charge in the active LAMMPS unit style.
    pub charge: Option<f64>,
    /// Coordinate in the active LAMMPS unit style.
    pub position: [f64; 3],
    /// Optional periodic image counters.
    pub image: Option<[i32; 3]>,
}

/// One typed bonded interaction retaining file atom identifiers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LammpsInteraction<const N: usize> {
    /// File-level interaction identifier.
    pub id: i64,
    /// Force-field interaction type.
    pub interaction_type: u32,
    /// Atom identifiers in declared order.
    pub atoms: [i64; N],
}

/// Triclinic-capable simulation cell from a data header.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LammpsDataCell {
    /// Lower bounds before tilt correction.
    pub lower: [f64; 3],
    /// Upper bounds before tilt correction.
    pub upper: [f64; 3],
    /// `xy`, `xz`, and `yz` tilt factors.
    pub tilt: [f64; 3],
}

/// Parsed LAMMPS data topology and force-field annotations.
#[derive(Clone, Debug, PartialEq)]
pub struct LammpsData {
    /// First line, retained verbatim.
    pub title: Box<str>,
    /// Atom style controlling record interpretation.
    pub atom_style: LammpsAtomStyle,
    /// Simulation cell when all three bounds are present.
    pub cell: Option<LammpsDataCell>,
    /// Mass by atom type.
    pub masses: BTreeMap<u32, f64>,
    /// Atoms sorted by file identifier.
    pub atoms: Vec<LammpsDataAtom>,
    /// Bond records sorted by interaction identifier.
    pub bonds: Vec<LammpsInteraction<2>>,
    /// Angle records sorted by interaction identifier.
    pub angles: Vec<LammpsInteraction<3>>,
    /// Dihedral records sorted by interaction identifier.
    pub dihedrals: Vec<LammpsInteraction<4>>,
    /// Improper records sorted by interaction identifier.
    pub impropers: Vec<LammpsInteraction<4>>,
    /// Unsupported sections retained as comment-free records.
    pub other_sections: BTreeMap<Box<str>, Vec<Box<str>>>,
}

/// Malformed or ambiguous LAMMPS data input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum LammpsDataError {
    /// The `Atoms` section was absent.
    #[error("LAMMPS data has no Atoms section")]
    MissingAtoms,
    /// The atom style was absent or unsupported.
    #[error("LAMMPS Atoms section must declare atomic, charge, molecular or full style")]
    MissingAtomStyle,
    /// A numeric header or section field was malformed.
    #[error("invalid LAMMPS data numeric field")]
    InvalidNumber,
    /// A record had the wrong number of fields for its section.
    #[error("invalid LAMMPS data record width")]
    InvalidRecord,
    /// A file-level identifier was repeated.
    #[error("duplicate LAMMPS data identifier")]
    DuplicateIdentifier,
    /// A connectivity record names an absent atom.
    #[error("LAMMPS connectivity references an absent atom")]
    UnknownAtom,
    /// A declared header count disagrees with the parsed section.
    #[error("LAMMPS data header count disagrees with its section")]
    CountMismatch,
}

/// Parses LAMMPS data files using `atomic`, `charge`, `molecular` or `full` atoms.
///
/// Coefficient and extension sections not interpreted here are retained in
/// [`LammpsData::other_sections`] rather than discarded.
///
/// # Errors
///
/// Returns an explicit error for ambiguous atom style, malformed numeric data,
/// duplicate identifiers, invalid connectivity, or count disagreement.
pub fn parse_lammps_data(text: &str) -> Result<LammpsData, LammpsDataError> {
    let mut lines = text.lines();
    let title: Box<str> = match lines.next() {
        Some(line) => line.into(),
        None => return Err(LammpsDataError::MissingAtoms),
    };
    let mut header = Header::default();
    let mut sections = BTreeMap::<Box<str>, Section>::new();
    let mut current: Option<Box<str>> = None;
    for raw in lines {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(section) = section_header(line) {
            current = Some(section.name.clone());
            sections.insert(section.name.clone(), section);
            continue;
        }
        let content = strip_comment(line);
        if content.is_empty() {
            continue;
        }
        match &current {
            Some(name) => {
                let Some(section) = sections.get_mut(name) else {
                    return Err(LammpsDataError::InvalidRecord);
                };
                section.records.push(content.into());
            }
            None => header.read(content)?,
        }
    }
    build(title, &header, sections)
}

#[derive(Default)]
struct Header {
    counts: BTreeMap<Box<str>, usize>,
    lower: [Option<f64>; 3],
    upper: [Option<f64>; 3],
    tilt: [f64; 3],
}

impl Header {
    fn read(&mut self, line: &str) -> Result<(), LammpsDataError> {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() >= 2
            && let Ok(count) = fields[0].parse::<usize>()
            && matches!(
                fields[1],
                "atoms" | "bonds" | "angles" | "dihedrals" | "impropers"
            )
        {
            self.counts.insert(fields[1].into(), count);
            return Ok(());
        }
        if fields.len() == 4 && fields[2].ends_with("lo") && fields[3].ends_with("hi") {
            let axis = match fields[2].as_bytes().first() {
                Some(b'x') => 0,
                Some(b'y') => 1,
                Some(b'z') => 2,
                _ => return Err(LammpsDataError::InvalidRecord),
            };
            self.lower[axis] = Some(number(fields[0])?);
            self.upper[axis] = Some(number(fields[1])?);
        } else if fields.len() == 6 && fields[3..] == ["xy", "xz", "yz"] {
            self.tilt = [number(fields[0])?, number(fields[1])?, number(fields[2])?];
        }
        Ok(())
    }
}

struct Section {
    name: Box<str>,
    style: Option<Box<str>>,
    records: Vec<Box<str>>,
}

fn section_header(line: &str) -> Option<Section> {
    let (name, style) = match line.split_once('#') {
        Some((name, style)) => (name.trim(), Some(style.trim().into())),
        None => (line, None),
    };
    let first = name.split_whitespace().next()?;
    if !first.as_bytes().first()?.is_ascii_uppercase() {
        return None;
    }
    Some(Section {
        name: name.into(),
        style,
        records: Vec::new(),
    })
}

fn strip_comment(line: &str) -> &str {
    match line.split_once('#') {
        Some((content, _)) => content.trim(),
        None => line,
    }
}

fn build(
    title: Box<str>,
    header: &Header,
    mut sections: BTreeMap<Box<str>, Section>,
) -> Result<LammpsData, LammpsDataError> {
    let atoms_section = sections
        .remove("Atoms")
        .ok_or(LammpsDataError::MissingAtoms)?;
    let atom_style = parse_style(atoms_section.style.as_deref())?;
    let mut atoms = parse_atoms(&atoms_section.records, atom_style)?;
    atoms.sort_unstable_by_key(|atom| atom.id);
    unique(atoms.iter().map(|atom| atom.id))?;
    let atom_ids: BTreeSet<_> = atoms.iter().map(|atom| atom.id).collect();
    let masses = parse_masses(sections.remove("Masses"))?;
    let mut bonds = parse_interactions(sections.remove("Bonds"), &atom_ids)?;
    let mut angles = parse_interactions(sections.remove("Angles"), &atom_ids)?;
    let mut dihedrals = parse_interactions(sections.remove("Dihedrals"), &atom_ids)?;
    let mut impropers = parse_interactions(sections.remove("Impropers"), &atom_ids)?;
    bonds.sort_unstable_by_key(|item| item.id);
    angles.sort_unstable_by_key(|item| item.id);
    dihedrals.sort_unstable_by_key(|item| item.id);
    impropers.sort_unstable_by_key(|item| item.id);
    validate_counts(
        &header.counts,
        atoms.len(),
        bonds.len(),
        angles.len(),
        dihedrals.len(),
        impropers.len(),
    )?;
    let cell = match (header.lower, header.upper) {
        ([Some(xl), Some(yl), Some(zl)], [Some(xu), Some(yu), Some(zu)]) => Some(LammpsDataCell {
            lower: [xl, yl, zl],
            upper: [xu, yu, zu],
            tilt: header.tilt,
        }),
        _ => None,
    };
    let other_sections = sections
        .into_iter()
        .map(|(name, section)| (name, section.records))
        .collect();
    Ok(LammpsData {
        title,
        atom_style,
        cell,
        masses,
        atoms,
        bonds,
        angles,
        dihedrals,
        impropers,
        other_sections,
    })
}

fn parse_style(style: Option<&str>) -> Result<LammpsAtomStyle, LammpsDataError> {
    match style {
        Some("atomic") => Ok(LammpsAtomStyle::Atomic),
        Some("charge") => Ok(LammpsAtomStyle::Charge),
        Some("molecular") => Ok(LammpsAtomStyle::Molecular),
        Some("full") => Ok(LammpsAtomStyle::Full),
        _ => Err(LammpsDataError::MissingAtomStyle),
    }
}

fn parse_atoms(
    records: &[Box<str>],
    style: LammpsAtomStyle,
) -> Result<Vec<LammpsDataAtom>, LammpsDataError> {
    records.iter().map(|line| parse_atom(line, style)).collect()
}

fn parse_atom(line: &str, style: LammpsAtomStyle) -> Result<LammpsDataAtom, LammpsDataError> {
    let fields: Vec<_> = line.split_whitespace().collect();
    let base = match style {
        LammpsAtomStyle::Atomic => 5,
        LammpsAtomStyle::Charge | LammpsAtomStyle::Molecular => 6,
        LammpsAtomStyle::Full => 7,
    };
    if fields.len() != base && fields.len() != base + 3 {
        return Err(LammpsDataError::InvalidRecord);
    }
    let id = integer(fields[0])?;
    let (molecule, atom_type_index, charge, position_index) = match style {
        LammpsAtomStyle::Atomic => (None, 1, None, 2),
        LammpsAtomStyle::Charge => (None, 1, Some(number(fields[2])?), 3),
        LammpsAtomStyle::Molecular => (Some(integer(fields[1])?), 2, None, 3),
        LammpsAtomStyle::Full => (Some(integer(fields[1])?), 2, Some(number(fields[3])?), 4),
    };
    let image = if fields.len() == base + 3 {
        Some([
            signed(fields[base])?,
            signed(fields[base + 1])?,
            signed(fields[base + 2])?,
        ])
    } else {
        None
    };
    Ok(LammpsDataAtom {
        id,
        molecule,
        atom_type: unsigned(fields[atom_type_index])?,
        charge,
        position: [
            number(fields[position_index])?,
            number(fields[position_index + 1])?,
            number(fields[position_index + 2])?,
        ],
        image,
    })
}

fn parse_masses(section: Option<Section>) -> Result<BTreeMap<u32, f64>, LammpsDataError> {
    let mut masses = BTreeMap::new();
    if let Some(section) = section {
        for line in section.records {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() < 2
                || masses
                    .insert(unsigned(fields[0])?, number(fields[1])?)
                    .is_some()
            {
                return Err(LammpsDataError::DuplicateIdentifier);
            }
        }
    }
    Ok(masses)
}

fn parse_interactions<const N: usize>(
    section: Option<Section>,
    atom_ids: &BTreeSet<i64>,
) -> Result<Vec<LammpsInteraction<N>>, LammpsDataError> {
    let Some(section) = section else {
        return Ok(Vec::new());
    };
    let mut output = Vec::with_capacity(section.records.len());
    for line in section.records {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != N + 2 {
            return Err(LammpsDataError::InvalidRecord);
        }
        let mut atoms = [0; N];
        for (index, atom) in atoms.iter_mut().enumerate() {
            *atom = integer(fields[index + 2])?;
        }
        if atoms.iter().any(|atom| !atom_ids.contains(atom)) {
            return Err(LammpsDataError::UnknownAtom);
        }
        output.push(LammpsInteraction {
            id: integer(fields[0])?,
            interaction_type: unsigned(fields[1])?,
            atoms,
        });
    }
    unique(output.iter().map(|item| item.id))?;
    Ok(output)
}

fn unique(values: impl IntoIterator<Item = i64>) -> Result<(), LammpsDataError> {
    let mut seen = BTreeSet::new();
    for value in values {
        if !seen.insert(value) {
            return Err(LammpsDataError::DuplicateIdentifier);
        }
    }
    Ok(())
}
fn validate_counts(
    counts: &BTreeMap<Box<str>, usize>,
    atoms: usize,
    bonds: usize,
    angles: usize,
    dihedrals: usize,
    impropers: usize,
) -> Result<(), LammpsDataError> {
    for (name, found) in [
        ("atoms", atoms),
        ("bonds", bonds),
        ("angles", angles),
        ("dihedrals", dihedrals),
        ("impropers", impropers),
    ] {
        if counts.get(name).is_some_and(|expected| *expected != found) {
            return Err(LammpsDataError::CountMismatch);
        }
    }
    Ok(())
}
fn number(value: &str) -> Result<f64, LammpsDataError> {
    value.parse().map_err(|_| LammpsDataError::InvalidNumber)
}
fn integer(value: &str) -> Result<i64, LammpsDataError> {
    value.parse().map_err(|_| LammpsDataError::InvalidNumber)
}
fn unsigned(value: &str) -> Result<u32, LammpsDataError> {
    value.parse().map_err(|_| LammpsDataError::InvalidNumber)
}
fn signed(value: &str) -> Result<i32, LammpsDataError> {
    value.parse().map_err(|_| LammpsDataError::InvalidNumber)
}

#[cfg(test)]
#[path = "lammps_data_tests.rs"]
mod tests;
