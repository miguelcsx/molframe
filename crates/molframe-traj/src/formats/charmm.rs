//! Faithful CHARMM coordinate CARD records, including `EXT` and `FREE`.

use crate::Frame;
use std::fmt::Write as _;

/// Coordinate-card field layout.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CharmmCardFormat {
    /// Original fixed-width fields.
    Standard,
    /// Expanded `EXT` fixed-width fields.
    Extended,
    /// Whitespace-delimited `FREE` fields.
    Free,
}

/// One complete CHARMM atom card.
#[derive(Clone, Debug, PartialEq)]
pub struct CharmmAtom {
    /// Atom serial.
    pub atom_number: i64,
    /// Residue serial.
    pub residue_number: i64,
    /// Residue name.
    pub residue_name: Box<str>,
    /// Atom name.
    pub atom_name: Box<str>,
    /// Cartesian position in ångström.
    pub position: [f32; 3],
    /// Segment identifier.
    pub segment_id: Box<str>,
    /// Residue identifier field.
    pub residue_id: Box<str>,
    /// Weighting scalar.
    pub weight: f32,
}

/// One coordinate CARD with title records and explicit layout.
#[derive(Clone, Debug, PartialEq)]
pub struct CharmmCard {
    /// Leading title records including their `*` prefix.
    pub titles: Vec<Box<str>>,
    /// Selected field layout.
    pub format: CharmmCardFormat,
    /// Atom cards in source order.
    pub atoms: Vec<CharmmAtom>,
}

impl CharmmCard {
    /// Returns the common coordinate-only frame.
    #[must_use]
    pub fn to_frame(&self) -> Frame {
        Frame {
            positions: self.atoms.iter().map(|atom| atom.position).collect(),
        }
    }
}

/// Invalid or unrepresentable coordinate CARD data.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CharmmError {
    /// Count record is absent or malformed.
    #[error("CHARMM CARD atom count is absent or invalid")]
    InvalidAtomCount,
    /// A fixed-width atom line is truncated.
    #[error("CHARMM atom card is too short")]
    ShortAtomCard,
    /// A required numeric or textual field is invalid.
    #[error("invalid CHARMM atom-card field")]
    InvalidField,
    /// Declared and parsed atom counts differ.
    #[error("CHARMM CARD expected {expected} atoms, found {found}")]
    AtomCountMismatch {
        /// Declared atom count.
        expected: usize,
        /// Parsed atom count.
        found: usize,
    },
    /// Selected layout cannot represent a value.
    #[error("CHARMM CARD value cannot be represented by its selected layout")]
    Unrepresentable,
}

/// Parses a complete coordinate CARD and retains every semantic atom field.
///
/// # Errors
///
/// Returns [`CharmmError`] for malformed counts, fields, or truncated atoms.
pub fn parse_charmm_record(text: &str) -> Result<CharmmCard, CharmmError> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let mut titles = Vec::new();
    let count_line = loop {
        let line = lines.next().ok_or(CharmmError::InvalidAtomCount)?;
        if line.trim_start().starts_with('*') {
            titles.push(line.into());
        } else {
            break line;
        }
    };
    let mut count_fields = count_line.split_whitespace();
    let atom_count = number(count_fields.next()).map_err(|_| CharmmError::InvalidAtomCount)?;
    let markers: Vec<_> = count_fields.collect();
    let format = if markers
        .iter()
        .any(|value| value.eq_ignore_ascii_case("FREE"))
    {
        CharmmCardFormat::Free
    } else if markers
        .iter()
        .any(|value| value.eq_ignore_ascii_case("EXT"))
    {
        CharmmCardFormat::Extended
    } else {
        CharmmCardFormat::Standard
    };
    let atoms = lines
        .take(atom_count)
        .map(|line| parse_atom(line, format))
        .collect::<Result<Vec<_>, _>>()?;
    if atoms.len() != atom_count {
        return Err(CharmmError::AtomCountMismatch {
            expected: atom_count,
            found: atoms.len(),
        });
    }
    Ok(CharmmCard {
        titles,
        format,
        atoms,
    })
}

/// Writes a complete coordinate CARD without synthesizing topology fields.
///
/// # Errors
///
/// Returns [`CharmmError::Unrepresentable`] when a selected width would overflow.
pub fn write_charmm_card(card: &CharmmCard) -> Result<String, CharmmError> {
    validate(card)?;
    let mut output = String::new();
    for title in &card.titles {
        writeln!(output, "{title}").map_err(|_| CharmmError::Unrepresentable)?;
    }
    match card.format {
        CharmmCardFormat::Standard => writeln!(output, "{:>5}", card.atoms.len()),
        CharmmCardFormat::Extended => writeln!(output, "{:>10} EXT", card.atoms.len()),
        CharmmCardFormat::Free => writeln!(output, "{} FREE", card.atoms.len()),
    }
    .map_err(|_| CharmmError::Unrepresentable)?;
    for atom in &card.atoms {
        write_atom(&mut output, atom, card.format)?;
    }
    Ok(output)
}

fn parse_atom(line: &str, format: CharmmCardFormat) -> Result<CharmmAtom, CharmmError> {
    match format {
        CharmmCardFormat::Free => parse_fields(&line.split_whitespace().collect::<Vec<_>>()),
        CharmmCardFormat::Standard => parse_fixed(
            line,
            [
                0, 5, 5, 10, 11, 15, 16, 20, 20, 30, 40, 50, 51, 55, 56, 60, 60, 70,
            ],
        ),
        CharmmCardFormat::Extended => parse_fixed(
            line,
            [
                0, 10, 10, 20, 22, 30, 32, 40, 40, 60, 80, 100, 102, 110, 112, 120, 120, 140,
            ],
        ),
    }
}

fn parse_fixed(line: &str, offsets: [usize; 18]) -> Result<CharmmAtom, CharmmError> {
    if line.len() < offsets[17] {
        return Err(CharmmError::ShortAtomCard);
    }
    let fields = [
        field(line, offsets[0], offsets[1])?,
        field(line, offsets[2], offsets[3])?,
        field(line, offsets[4], offsets[5])?,
        field(line, offsets[6], offsets[7])?,
        field(line, offsets[8], offsets[9])?,
        field(line, offsets[9], offsets[10])?,
        field(line, offsets[10], offsets[11])?,
        field(line, offsets[12], offsets[13])?,
        field(line, offsets[14], offsets[15])?,
        field(line, offsets[16], offsets[17])?,
    ];
    parse_fields(&fields)
}

fn parse_fields(fields: &[&str]) -> Result<CharmmAtom, CharmmError> {
    if fields.len() < 10 {
        return Err(CharmmError::ShortAtomCard);
    }
    Ok(CharmmAtom {
        atom_number: number(Some(fields[0]))?,
        residue_number: number(Some(fields[1]))?,
        residue_name: fields[2].trim().into(),
        atom_name: fields[3].trim().into(),
        position: [
            number(Some(fields[4]))?,
            number(Some(fields[5]))?,
            number(Some(fields[6]))?,
        ],
        segment_id: fields[7].trim().into(),
        residue_id: fields[8].trim().into(),
        weight: number(Some(fields[9]))?,
    })
}

fn write_atom(
    output: &mut String,
    atom: &CharmmAtom,
    format: CharmmCardFormat,
) -> Result<(), CharmmError> {
    match format {
        CharmmCardFormat::Standard => writeln!(
            output,
            "{:>5}{:>5} {:<4} {:<4}{:>10.5}{:>10.5}{:>10.5} {:<4} {:<4}{:>10.5}",
            atom.atom_number,
            atom.residue_number,
            atom.residue_name,
            atom.atom_name,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            atom.segment_id,
            atom.residue_id,
            atom.weight
        ),
        CharmmCardFormat::Extended => writeln!(
            output,
            "{:>10}{:>10}  {:<8}  {:<8}{:>20.10}{:>20.10}{:>20.10}  {:<8}  {:<8}{:>20.10}",
            atom.atom_number,
            atom.residue_number,
            atom.residue_name,
            atom.atom_name,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            atom.segment_id,
            atom.residue_id,
            atom.weight
        ),
        CharmmCardFormat::Free => writeln!(
            output,
            "{} {} {} {} {:.10} {:.10} {:.10} {} {} {:.10}",
            atom.atom_number,
            atom.residue_number,
            atom.residue_name,
            atom.atom_name,
            atom.position[0],
            atom.position[1],
            atom.position[2],
            atom.segment_id,
            atom.residue_id,
            atom.weight
        ),
    }
    .map_err(|_| CharmmError::Unrepresentable)
}

fn validate(card: &CharmmCard) -> Result<(), CharmmError> {
    let width = match card.format {
        CharmmCardFormat::Standard => 4,
        CharmmCardFormat::Extended => 8,
        CharmmCardFormat::Free => usize::MAX,
    };
    if card
        .titles
        .iter()
        .any(|title| title.contains(['\n', '\r']) || !title.trim_start().starts_with('*'))
        || card.atoms.iter().any(|atom| {
            [
                atom.residue_name.as_ref(),
                atom.atom_name.as_ref(),
                atom.segment_id.as_ref(),
                atom.residue_id.as_ref(),
            ]
            .iter()
            .any(|value| {
                value.is_empty() || value.len() > width || value.chars().any(char::is_whitespace)
            }) || atom.position.iter().any(|value| !value.is_finite())
                || !atom.weight.is_finite()
        })
    {
        Err(CharmmError::Unrepresentable)
    } else {
        Ok(())
    }
}

fn field(line: &str, start: usize, end: usize) -> Result<&str, CharmmError> {
    line.get(start..end)
        .map(str::trim)
        .ok_or(CharmmError::ShortAtomCard)
}
fn number<T: std::str::FromStr>(value: Option<&str>) -> Result<T, CharmmError> {
    value
        .and_then(|value| value.trim().parse().ok())
        .ok_or(CharmmError::InvalidField)
}

#[cfg(test)]
#[path = "charmm_tests.rs"]
mod tests;
