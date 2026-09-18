//! Normal-mode data files, the interchange format for precomputed modes.
//!
//! A mode set computed elsewhere — by a dedicated NMA package, from an
//! essential-dynamics analysis of a long trajectory, or from experiment — is
//! the same thing [`crate::AnisotropicNetworkModel`] produces, and a caller
//! should be able to use one wherever the other fits. This reads and writes
//! that interchange form.
//!
//! The format is line-oriented: a keyword, then whitespace-separated values to
//! the end of the line. Coordinates and every mode carry three values per site,
//! so a file's shape is validated against its own coordinate block rather than
//! trusted. Reading is one pass, `O(values)`.

use std::fmt::Write as _;

/// One named displacement field over the file's sites.
#[derive(Clone, Debug, PartialEq)]
pub struct NormalMode {
    /// Mode index as written, which is not always its position in the file.
    pub index: i32,
    /// Scale factor the writer attached, conventionally the inverse square
    /// root of the eigenvalue so that a mode's drawn length reflects how
    /// easily it is excited.
    pub scale: f64,
    /// One displacement per site.
    pub displacements: Vec<[f64; 3]>,
}

/// A parsed normal-mode data file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NormalModeSet {
    /// Title, when the file declared one.
    pub name: Option<Box<str>>,
    /// Atom names, one per site, when present.
    pub atom_names: Vec<Box<str>>,
    /// Residue names, one per site, when present.
    pub residue_names: Vec<Box<str>>,
    /// Residue numbers, one per site, when present.
    pub residue_ids: Vec<i32>,
    /// Chain identifiers, one per site, when present.
    pub chain_ids: Vec<Box<str>>,
    /// Temperature factors, one per site, when present.
    pub b_factors: Vec<f64>,
    /// Reference coordinates the modes displace.
    pub coordinates: Vec<[f64; 3]>,
    /// Modes in file order.
    pub modes: Vec<NormalMode>,
}

impl NormalModeSet {
    /// Sites the file describes, taken from its coordinate block.
    #[must_use]
    pub fn len(&self) -> usize {
        self.coordinates.len()
    }

    /// Whether the file described no sites.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.coordinates.is_empty()
    }

    /// Writes `coordinates + sum(amplitude x mode)` into `out`.
    ///
    /// Amplitudes are applied as given; multiplying by each mode's own `scale`
    /// first is the caller's choice, because whether a drawn mode should be
    /// weighted by its excitability depends on what the picture is for.
    /// `O(modes x sites)`.
    pub fn displace(&self, amplitudes: &[f64], out: &mut [[f64; 3]]) {
        let count = self.coordinates.len().min(out.len());
        out[..count].copy_from_slice(&self.coordinates[..count]);
        for (amplitude, mode) in amplitudes.iter().zip(&self.modes) {
            if *amplitude == 0.0 {
                continue;
            }
            for (position, displacement) in out[..count].iter_mut().zip(&mode.displacements) {
                for (axis, offset) in position.iter_mut().zip(displacement) {
                    *axis += amplitude * offset;
                }
            }
        }
    }
}

/// Why a normal-mode data file could not be read.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum NmdError {
    /// The text is not UTF-8.
    #[error("normal-mode data is not valid UTF-8")]
    NotUtf8,
    /// A record held a token that is not a number.
    #[error("normal-mode {keyword} record holds a value that is not a number")]
    InvalidNumber {
        /// The record's keyword.
        keyword: &'static str,
    },
    /// A vector record's length is not a multiple of three.
    #[error("normal-mode {keyword} record is not a whole number of vectors")]
    UnalignedVectors {
        /// The record's keyword.
        keyword: &'static str,
    },
    /// A mode covers a different number of sites than the coordinates.
    #[error("normal-mode {found}-site mode does not match the {expected}-site coordinates")]
    SiteMismatch {
        /// Sites the mode covered.
        found: usize,
        /// Sites the coordinate block established.
        expected: usize,
    },
    /// A mode appeared before the coordinates that give it meaning.
    #[error("a normal mode appeared before the coordinate block")]
    ModeBeforeCoordinates,
}

/// Reads a normal-mode data file.
///
/// # Errors
///
/// Returns [`NmdError`] for non-UTF-8 text, a non-numeric or misaligned
/// record, or a mode whose site count disagrees with the coordinates.
pub fn read_nmd(bytes: &[u8]) -> Result<NormalModeSet, NmdError> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return Err(NmdError::NotUtf8);
    };
    let mut set = NormalModeSet::default();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split_ascii_whitespace();
        let Some(keyword) = fields.next() else {
            continue;
        };
        match keyword {
            "name" | "title" => {
                let rest: Vec<&str> = fields.collect();
                if !rest.is_empty() {
                    set.name = Some(Box::from(rest.join(" ")));
                }
            }
            "atomnames" => set.atom_names = fields.map(Box::from).collect(),
            "resnames" => set.residue_names = fields.map(Box::from).collect(),
            "chainids" | "chids" => set.chain_ids = fields.map(Box::from).collect(),
            "resids" => set.residue_ids = parse_i32(fields, "resids")?,
            "bfactors" | "betas" => set.b_factors = parse_f64(fields, "bfactors")?,
            "coordinates" => {
                set.coordinates = parse_vectors(fields, "coordinates")?;
            }
            "mode" => {
                let mode = read_mode(fields, set.coordinates.len())?;
                set.modes.push(mode);
            }
            // Everything else is viewer styling the model does not depend on.
            _ => {}
        }
    }
    Ok(set)
}

/// Writes a normal-mode data file.
///
/// Records the writer has nothing to say about are omitted rather than written
/// empty, so a round trip through this pair is stable.
#[must_use]
pub fn write_nmd(set: &NormalModeSet) -> String {
    let mut text = String::new();
    if let Some(name) = &set.name {
        let _ = writeln!(text, "name {name}");
    }
    write_labels(&mut text, "atomnames", &set.atom_names);
    write_labels(&mut text, "resnames", &set.residue_names);
    write_labels(&mut text, "chainids", &set.chain_ids);
    if !set.residue_ids.is_empty() {
        let _ = write!(text, "resids");
        for id in &set.residue_ids {
            let _ = write!(text, " {id}");
        }
        text.push('\n');
    }
    if !set.b_factors.is_empty() {
        let _ = write!(text, "bfactors");
        for value in &set.b_factors {
            let _ = write!(text, " {value}");
        }
        text.push('\n');
    }
    if !set.coordinates.is_empty() {
        let _ = write!(text, "coordinates");
        for point in &set.coordinates {
            for value in point {
                let _ = write!(text, " {value}");
            }
        }
        text.push('\n');
    }
    for mode in &set.modes {
        let _ = write!(text, "mode {} {}", mode.index, mode.scale);
        for displacement in &mode.displacements {
            for value in displacement {
                let _ = write!(text, " {value}");
            }
        }
        text.push('\n');
    }
    text
}

fn write_labels(text: &mut String, keyword: &str, labels: &[Box<str>]) {
    if labels.is_empty() {
        return;
    }
    let _ = write!(text, "{keyword}");
    for label in labels {
        let _ = write!(text, " {label}");
    }
    text.push('\n');
}

/// A mode record opens with an index and a scale, then its displacements.
///
/// Some writers omit one or both leading numbers; a record whose value count is
/// already a whole number of vectors is taken to have done so.
fn read_mode<'a>(
    fields: impl Iterator<Item = &'a str>,
    expected: usize,
) -> Result<NormalMode, NmdError> {
    if expected == 0 {
        return Err(NmdError::ModeBeforeCoordinates);
    }
    let tokens: Vec<&'a str> = fields.collect();
    let leading = tokens.len().saturating_sub(expected * 3);
    if tokens.len() < expected * 3 {
        return Err(NmdError::SiteMismatch {
            found: tokens.len() / 3,
            expected,
        });
    }
    let mut index = 0;
    let mut scale = 1.0;
    if leading >= 1 {
        index = match tokens[0].parse::<i32>() {
            Ok(value) => value,
            Err(_) => return Err(NmdError::InvalidNumber { keyword: "mode" }),
        };
    }
    if leading >= 2 {
        scale = match tokens[1].parse::<f64>() {
            Ok(value) => value,
            Err(_) => return Err(NmdError::InvalidNumber { keyword: "mode" }),
        };
    }
    let displacements = parse_vectors(tokens[leading..].iter().copied(), "mode")?;
    if displacements.len() != expected {
        return Err(NmdError::SiteMismatch {
            found: displacements.len(),
            expected,
        });
    }
    Ok(NormalMode {
        index,
        scale,
        displacements,
    })
}

fn parse_vectors<'a>(
    fields: impl Iterator<Item = &'a str>,
    keyword: &'static str,
) -> Result<Vec<[f64; 3]>, NmdError> {
    let values = parse_f64(fields, keyword)?;
    if !values.len().is_multiple_of(3) {
        return Err(NmdError::UnalignedVectors { keyword });
    }
    Ok(values
        .as_chunks::<3>()
        .0
        .iter()
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect())
}

fn parse_f64<'a>(
    fields: impl Iterator<Item = &'a str>,
    keyword: &'static str,
) -> Result<Vec<f64>, NmdError> {
    let mut values = Vec::new();
    for token in fields {
        match token.parse::<f64>() {
            Ok(value) => values.push(value),
            Err(_) => return Err(NmdError::InvalidNumber { keyword }),
        }
    }
    Ok(values)
}

fn parse_i32<'a>(
    fields: impl Iterator<Item = &'a str>,
    keyword: &'static str,
) -> Result<Vec<i32>, NmdError> {
    let mut values = Vec::new();
    for token in fields {
        match token.parse::<i32>() {
            Ok(value) => values.push(value),
            Err(_) => return Err(NmdError::InvalidNumber { keyword }),
        }
    }
    Ok(values)
}

#[cfg(test)]
#[path = "nmd_tests.rs"]
mod tests;
