//! AMBER `%FLAG` parameter/topology reader.
//!
//! Fixed-width fields are decoded from their declared Fortran format rather
//! than split on whitespace. This preserves adjacent negative values and the
//! four-character atom/residue fields used by canonical `parm7` writers.

use std::collections::BTreeMap;

const AMBER_CHARGE_SCALE: f64 = 18.2223;
type SectionMap = BTreeMap<Box<str>, AmberSection>;

/// One raw, self-describing `%FLAG` section.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmberSection {
    /// Original `%FORMAT` body.
    pub format: Box<str>,
    /// Fixed-width values with padding removed.
    pub values: Vec<Box<str>>,
}

/// Atom annotations joined from the parallel AMBER topology arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct AmberTopologyAtom {
    /// Atom name.
    pub name: Box<str>,
    /// Charge in proton-charge units.
    pub charge: f64,
    /// Mass in daltons.
    pub mass: f64,
    /// One-based force-field type index.
    pub type_index: u32,
    /// AMBER atom-type label.
    pub atom_type: Box<str>,
    /// Atomic number when the section is present.
    pub atomic_number: Option<u8>,
    /// Zero-based residue index.
    pub residue: u32,
}

/// Residue name and half-open atom range.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmberTopologyResidue {
    /// Residue label.
    pub name: Box<str>,
    /// First zero-based atom index.
    pub atom_start: u32,
    /// Exclusive zero-based atom end.
    pub atom_end: u32,
}

/// Bond and its one-based force-field parameter index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AmberTopologyBond {
    /// Zero-based atom indices.
    pub atoms: [u32; 2],
    /// One-based parameter-table index.
    pub parameter: u32,
    /// Whether the source section included hydrogen.
    pub includes_hydrogen: bool,
}

/// Angle and its one-based force-field parameter index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AmberTopologyAngle {
    /// Zero-based atom indices.
    pub atoms: [u32; 3],
    /// One-based parameter-table index.
    pub parameter: u32,
    /// Whether the source section included hydrogen.
    pub includes_hydrogen: bool,
}

/// Dihedral/improper and its AMBER sign semantics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AmberTopologyDihedral {
    /// Zero-based atom indices after decoding coordinate-array offsets.
    pub atoms: [u32; 4],
    /// One-based parameter-table index.
    pub parameter: u32,
    /// Negative fourth encoded atom marks an improper.
    pub improper: bool,
    /// Negative third encoded atom suppresses the 1-4 nonbonded interaction.
    pub ignore_end_group: bool,
    /// Whether the source section included hydrogen.
    pub includes_hydrogen: bool,
}

/// Fully joined AMBER topology plus every raw section.
#[derive(Clone, Debug, PartialEq)]
pub struct AmberTopology {
    /// `%VERSION` line when present.
    pub version: Option<Box<str>>,
    /// Joined atom table.
    pub atoms: Vec<AmberTopologyAtom>,
    /// Joined residue table.
    pub residues: Vec<AmberTopologyResidue>,
    /// Bonds from both hydrogen partitions.
    pub bonds: Vec<AmberTopologyBond>,
    /// Angles from both hydrogen partitions.
    pub angles: Vec<AmberTopologyAngle>,
    /// Proper and improper dihedrals from both partitions.
    pub dihedrals: Vec<AmberTopologyDihedral>,
    /// All sections, including force constants and unknown extensions.
    pub sections: BTreeMap<Box<str>, AmberSection>,
}

/// Malformed or inconsistent AMBER topology.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AmberTopologyError {
    /// A `%FLAG` did not have a following `%FORMAT`.
    #[error("AMBER topology flag has no valid format")]
    MissingFormat,
    /// A Fortran format was unsupported or malformed.
    #[error("invalid AMBER topology Fortran format")]
    InvalidFormat,
    /// A required parallel array was absent.
    #[error("AMBER topology is missing a required section")]
    MissingSection,
    /// A numeric field could not be parsed.
    #[error("invalid AMBER topology numeric field")]
    InvalidNumber,
    /// Parallel arrays or a pointer count disagree.
    #[error("AMBER topology section lengths are inconsistent")]
    LengthMismatch,
    /// A residue pointer or bonded atom offset is invalid.
    #[error("AMBER topology contains an invalid atom reference")]
    InvalidAtomReference,
}

/// Parses AMBER version-7 `prmtop`, `parm7`, and `%FLAG`-style `.top` files.
///
/// Charges are converted from AMBER's scaled representation to electron charge.
/// Every section remains available verbatim-at-field-level for force-field data
/// that the joined topology does not interpret.
///
/// # Errors
///
/// Returns an explicit error for malformed format declarations, missing or
/// inconsistent required arrays, invalid numbers, and bad atom references.
pub fn parse_amber_topology(text: &str) -> Result<AmberTopology, AmberTopologyError> {
    let (version, sections) = parse_sections(text)?;
    let names = required(&sections, "ATOM_NAME")?;
    let atom_count = names.values.len();
    validate_pointer_count(&sections, atom_count)?;
    let charges = floats(required(&sections, "CHARGE")?)?;
    let masses = floats(required(&sections, "MASS")?)?;
    let type_indices = unsigneds(required(&sections, "ATOM_TYPE_INDEX")?)?;
    let atom_types = required(&sections, "AMBER_ATOM_TYPE")?;
    let atomic_numbers = optional_unsigneds(&sections, "ATOMIC_NUMBER")?;
    for length in [
        charges.len(),
        masses.len(),
        type_indices.len(),
        atom_types.values.len(),
    ] {
        if length != atom_count {
            return Err(AmberTopologyError::LengthMismatch);
        }
    }
    if atomic_numbers
        .as_ref()
        .is_some_and(|values| values.len() != atom_count)
    {
        return Err(AmberTopologyError::LengthMismatch);
    }
    let residues = residues(&sections, atom_count)?;
    let residue_by_atom = residue_membership(&residues, atom_count)?;
    let atoms = (0..atom_count)
        .map(|index| AmberTopologyAtom {
            name: names.values[index].clone(),
            charge: charges[index] / AMBER_CHARGE_SCALE,
            mass: masses[index],
            type_index: type_indices[index],
            atom_type: atom_types.values[index].clone(),
            atomic_number: atomic_numbers
                .as_ref()
                .and_then(|values| values[index].try_into().ok()),
            residue: residue_by_atom[index],
        })
        .collect();
    let mut bonds = parse_bonds(&sections, "BONDS_INC_HYDROGEN", true, atom_count)?;
    bonds.extend(parse_bonds(
        &sections,
        "BONDS_WITHOUT_HYDROGEN",
        false,
        atom_count,
    )?);
    let mut angles = parse_angles(&sections, "ANGLES_INC_HYDROGEN", true, atom_count)?;
    angles.extend(parse_angles(
        &sections,
        "ANGLES_WITHOUT_HYDROGEN",
        false,
        atom_count,
    )?);
    let mut dihedrals = parse_dihedrals(&sections, "DIHEDRALS_INC_HYDROGEN", true, atom_count)?;
    dihedrals.extend(parse_dihedrals(
        &sections,
        "DIHEDRALS_WITHOUT_HYDROGEN",
        false,
        atom_count,
    )?);
    Ok(AmberTopology {
        version,
        atoms,
        residues,
        bonds,
        angles,
        dihedrals,
        sections,
    })
}

fn parse_sections(text: &str) -> Result<(Option<Box<str>>, SectionMap), AmberTopologyError> {
    let lines: Vec<_> = text.lines().collect();
    let version = lines
        .first()
        .filter(|line| line.starts_with("%VERSION"))
        .map(|line| (*line).into());
    let mut sections = BTreeMap::new();
    let mut cursor = usize::from(version.is_some());
    while cursor < lines.len() {
        if !lines[cursor].starts_with("%FLAG ") {
            cursor += 1;
            continue;
        }
        let name: Box<str> = lines[cursor][6..].trim().into();
        cursor += 1;
        let format_line = lines.get(cursor).ok_or(AmberTopologyError::MissingFormat)?;
        let format = format_line
            .trim()
            .strip_prefix("%FORMAT(")
            .and_then(|value| value.strip_suffix(')'))
            .ok_or(AmberTopologyError::MissingFormat)?;
        let layout = Layout::parse(format)?;
        cursor += 1;
        let mut values = Vec::new();
        while cursor < lines.len() && !lines[cursor].starts_with("%FLAG ") {
            layout.read_line(lines[cursor], &mut values);
            cursor += 1;
        }
        sections.insert(
            name,
            AmberSection {
                format: format.into(),
                values,
            },
        );
    }
    Ok((version, sections))
}

struct Layout {
    repeat: usize,
    width: usize,
}

impl Layout {
    fn parse(format: &str) -> Result<Self, AmberTopologyError> {
        let compact: String = format
            .chars()
            .filter(|character| !character.is_whitespace())
            .collect();
        let type_at = compact
            .find(|character: char| character.is_ascii_alphabetic())
            .ok_or(AmberTopologyError::InvalidFormat)?;
        let repeat = compact[..type_at]
            .parse()
            .map_err(|_| AmberTopologyError::InvalidFormat)?;
        let rest = &compact[type_at + 1..];
        let width_text = match rest.split_once('.') {
            Some((width, _)) => width,
            None => rest,
        };
        let width = width_text
            .parse()
            .map_err(|_| AmberTopologyError::InvalidFormat)?;
        if repeat == 0 || width == 0 {
            return Err(AmberTopologyError::InvalidFormat);
        }
        Ok(Self { repeat, width })
    }

    fn read_line(&self, line: &str, output: &mut Vec<Box<str>>) {
        for field in line.as_bytes().chunks(self.width).take(self.repeat) {
            let value = String::from_utf8_lossy(field);
            let value = value.trim();
            if !value.is_empty() {
                output.push(value.into());
            }
        }
    }
}

fn required<'a>(
    sections: &'a BTreeMap<Box<str>, AmberSection>,
    name: &str,
) -> Result<&'a AmberSection, AmberTopologyError> {
    sections.get(name).ok_or(AmberTopologyError::MissingSection)
}
fn floats(section: &AmberSection) -> Result<Vec<f64>, AmberTopologyError> {
    section
        .values
        .iter()
        .map(|value| {
            value
                .replace(['D', 'd'], "E")
                .parse()
                .map_err(|_| AmberTopologyError::InvalidNumber)
        })
        .collect()
}
fn unsigneds(section: &AmberSection) -> Result<Vec<u32>, AmberTopologyError> {
    section
        .values
        .iter()
        .map(|value| value.parse().map_err(|_| AmberTopologyError::InvalidNumber))
        .collect()
}
fn signeds(section: &AmberSection) -> Result<Vec<i64>, AmberTopologyError> {
    section
        .values
        .iter()
        .map(|value| value.parse().map_err(|_| AmberTopologyError::InvalidNumber))
        .collect()
}
fn optional_unsigneds(
    sections: &BTreeMap<Box<str>, AmberSection>,
    name: &str,
) -> Result<Option<Vec<u32>>, AmberTopologyError> {
    sections.get(name).map(unsigneds).transpose()
}

fn validate_pointer_count(
    sections: &BTreeMap<Box<str>, AmberSection>,
    atom_count: usize,
) -> Result<(), AmberTopologyError> {
    let pointers = unsigneds(required(sections, "POINTERS")?)?;
    if pointers.first().copied() != u32::try_from(atom_count).ok() {
        return Err(AmberTopologyError::LengthMismatch);
    }
    Ok(())
}

fn residues(
    sections: &BTreeMap<Box<str>, AmberSection>,
    atom_count: usize,
) -> Result<Vec<AmberTopologyResidue>, AmberTopologyError> {
    let labels = required(sections, "RESIDUE_LABEL")?;
    let pointers = unsigneds(required(sections, "RESIDUE_POINTER")?)?;
    if labels.values.len() != pointers.len() {
        return Err(AmberTopologyError::LengthMismatch);
    }
    let mut output = Vec::with_capacity(pointers.len());
    for index in 0..pointers.len() {
        let start = pointers[index]
            .checked_sub(1)
            .ok_or(AmberTopologyError::InvalidAtomReference)?;
        let end = if index + 1 < pointers.len() {
            pointers[index + 1]
                .checked_sub(1)
                .ok_or(AmberTopologyError::InvalidAtomReference)?
        } else {
            u32::try_from(atom_count).map_err(|_| AmberTopologyError::LengthMismatch)?
        };
        if start > end || usize::try_from(end).map_or(true, |end_index| end_index > atom_count) {
            return Err(AmberTopologyError::InvalidAtomReference);
        }
        output.push(AmberTopologyResidue {
            name: labels.values[index].clone(),
            atom_start: start,
            atom_end: end,
        });
    }
    Ok(output)
}

fn residue_membership(
    residues: &[AmberTopologyResidue],
    atom_count: usize,
) -> Result<Vec<u32>, AmberTopologyError> {
    let mut membership = vec![u32::MAX; atom_count];
    for (index, residue) in residues.iter().enumerate() {
        let residue_index = u32::try_from(index).map_err(|_| AmberTopologyError::LengthMismatch)?;
        for atom in residue.atom_start..residue.atom_end {
            let atom =
                usize::try_from(atom).map_err(|_| AmberTopologyError::InvalidAtomReference)?;
            membership[atom] = residue_index;
        }
    }
    if membership.contains(&u32::MAX) {
        return Err(AmberTopologyError::InvalidAtomReference);
    }
    Ok(membership)
}

fn encoded_atom(value: i64, atom_count: usize) -> Result<u32, AmberTopologyError> {
    let offset = value.unsigned_abs();
    let atom_count = u64::try_from(atom_count).map_err(|_| AmberTopologyError::LengthMismatch)?;
    if !offset.is_multiple_of(3) || offset / 3 >= atom_count {
        return Err(AmberTopologyError::InvalidAtomReference);
    }
    u32::try_from(offset / 3).map_err(|_| AmberTopologyError::InvalidAtomReference)
}

fn parse_bonds(
    sections: &BTreeMap<Box<str>, AmberSection>,
    name: &str,
    includes_hydrogen: bool,
    atom_count: usize,
) -> Result<Vec<AmberTopologyBond>, AmberTopologyError> {
    let Some(section) = sections.get(name) else {
        return Ok(Vec::new());
    };
    let values = signeds(section)?;
    if values.len() % 3 != 0 {
        return Err(AmberTopologyError::LengthMismatch);
    }
    values
        .chunks_exact(3)
        .map(|row| {
            Ok(AmberTopologyBond {
                atoms: [
                    encoded_atom(row[0], atom_count)?,
                    encoded_atom(row[1], atom_count)?,
                ],
                parameter: u32::try_from(row[2]).map_err(|_| AmberTopologyError::InvalidNumber)?,
                includes_hydrogen,
            })
        })
        .collect()
}

fn parse_angles(
    sections: &BTreeMap<Box<str>, AmberSection>,
    name: &str,
    includes_hydrogen: bool,
    atom_count: usize,
) -> Result<Vec<AmberTopologyAngle>, AmberTopologyError> {
    let Some(section) = sections.get(name) else {
        return Ok(Vec::new());
    };
    let values = signeds(section)?;
    if values.len() % 4 != 0 {
        return Err(AmberTopologyError::LengthMismatch);
    }
    values
        .chunks_exact(4)
        .map(|row| {
            Ok(AmberTopologyAngle {
                atoms: [
                    encoded_atom(row[0], atom_count)?,
                    encoded_atom(row[1], atom_count)?,
                    encoded_atom(row[2], atom_count)?,
                ],
                parameter: u32::try_from(row[3]).map_err(|_| AmberTopologyError::InvalidNumber)?,
                includes_hydrogen,
            })
        })
        .collect()
}

fn parse_dihedrals(
    sections: &BTreeMap<Box<str>, AmberSection>,
    name: &str,
    includes_hydrogen: bool,
    atom_count: usize,
) -> Result<Vec<AmberTopologyDihedral>, AmberTopologyError> {
    let Some(section) = sections.get(name) else {
        return Ok(Vec::new());
    };
    let values = signeds(section)?;
    if values.len() % 5 != 0 {
        return Err(AmberTopologyError::LengthMismatch);
    }
    values
        .chunks_exact(5)
        .map(|row| {
            Ok(AmberTopologyDihedral {
                atoms: [
                    encoded_atom(row[0], atom_count)?,
                    encoded_atom(row[1], atom_count)?,
                    encoded_atom(row[2], atom_count)?,
                    encoded_atom(row[3], atom_count)?,
                ],
                parameter: u32::try_from(row[4]).map_err(|_| AmberTopologyError::InvalidNumber)?,
                improper: row[3] < 0,
                ignore_end_group: row[2] < 0,
                includes_hydrogen,
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "amber_topology_tests.rs"]
mod tests;
