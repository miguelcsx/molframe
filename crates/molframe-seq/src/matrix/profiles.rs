//! Named score-table profiles parsed from the bundled NCBI BLAST archive.

use super::SubstitutionMatrix;

const BUNDLE: &str = include_str!("../../data/ncbi_matrices.txt");
const VERSION: &str = "ncbi-blast-1997-08-26";
const RETRIEVED: &str = "2026-08-08";
const SOURCE_ROOT: &str = "https://ftp.ncbi.nih.gov/blast/matrices/";

/// Stable provenance of one substitution score table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatrixIdentity {
    name: Box<str>,
    version: Box<str>,
    source: Box<str>,
    retrieved: Box<str>,
}

impl MatrixIdentity {
    pub(super) fn custom() -> Self {
        Self {
            name: "custom".into(),
            version: "caller-supplied".into(),
            source: "caller".into(),
            retrieved: "not-applicable".into(),
        }
    }

    /// NCBI profile name or `custom`.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Exact archive version used by the bundled profile.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Direct source URL for the table.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Date on which the bundled source was retrieved.
    #[must_use]
    pub fn retrieved(&self) -> &str {
        &self.retrieved
    }
}

/// Named substitution matrix family and level.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MatrixProfile {
    /// BLOSUM clustering level from 30 through 90 in increments of five.
    Blosum(u8),
    /// PAM evolutionary distance from 10 through 500 in increments of ten.
    Pam(u16),
    /// NCBI's exact protein identity table.
    Identity,
    /// NCBI `NUC.4.4` over the IUPAC nucleotide alphabet.
    Nuc44,
}

impl MatrixProfile {
    fn name(self) -> Result<String, MatrixError> {
        match self {
            Self::Blosum(level) if (30..=90).contains(&level) && level % 5 == 0 => {
                Ok(format!("BLOSUM{level}"))
            }
            Self::Pam(distance) if (10..=500).contains(&distance) && distance % 10 == 0 => {
                Ok(format!("PAM{distance}"))
            }
            Self::Identity => Ok("IDENTITY".to_owned()),
            Self::Nuc44 => Ok("NUC.4.4".to_owned()),
            Self::Blosum(_) | Self::Pam(_) => Err(MatrixError::UnknownProfile),
        }
    }
}

/// A named profile is unavailable or its bundled table is malformed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MatrixError {
    /// Family parameter is outside the bundled, versioned range.
    UnknownProfile,
    /// A bundled section is absent or does not form a square integer matrix.
    MalformedProfile,
    /// A table has no explicit ambiguity symbol for out-of-alphabet bytes.
    MissingUnknownSymbol,
}

impl std::fmt::Display for MatrixError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownProfile => formatter.write_str("substitution matrix profile unavailable"),
            Self::MalformedProfile => formatter.write_str("bundled substitution matrix malformed"),
            Self::MissingUnknownSymbol => {
                formatter.write_str("substitution matrix has no X or N ambiguity symbol")
            }
        }
    }
}

impl std::error::Error for MatrixError {}

/// Loads an inspectable, versioned NCBI substitution matrix.
///
/// # Errors
///
/// Returns [`MatrixError`] for an unsupported family parameter or invalid
/// bundled source data.
pub fn load_matrix(profile: MatrixProfile) -> Result<SubstitutionMatrix, MatrixError> {
    let name = profile.name()?;
    let section = section(&name).ok_or(MatrixError::MalformedProfile)?;
    let (alphabet, scores) = parse(section)?;
    let identity = MatrixIdentity {
        name: name.clone().into(),
        version: VERSION.into(),
        source: format!("{SOURCE_ROOT}{name}").into(),
        retrieved: RETRIEVED.into(),
    };
    SubstitutionMatrix::from_profile(&alphabet, scores, identity)
}

fn section(name: &str) -> Option<&'static str> {
    let marker = format!("@@ {name}\n");
    let start = BUNDLE.find(&marker)? + marker.len();
    let remaining = &BUNDLE[start..];
    let end = match remaining.find("\n@@ ") {
        Some(end) => end,
        None => remaining.len(),
    };
    remaining.get(..end)
}

fn parse(text: &str) -> Result<(Vec<u8>, Vec<i32>), MatrixError> {
    let mut rows = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'));
    let header = rows.next().ok_or(MatrixError::MalformedProfile)?;
    let alphabet = header
        .split_ascii_whitespace()
        .map(|symbol| {
            let bytes = symbol.as_bytes();
            match bytes {
                [byte] => Ok(*byte),
                _ => Err(MatrixError::MalformedProfile),
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut scores = Vec::with_capacity(alphabet.len() * alphabet.len());
    for (row_index, line) in rows.enumerate() {
        let mut fields = line.split_ascii_whitespace();
        let row_name = fields.next().ok_or(MatrixError::MalformedProfile)?;
        let Some(&expected) = alphabet.get(row_index) else {
            return Err(MatrixError::MalformedProfile);
        };
        if row_name.as_bytes() != [expected] {
            return Err(MatrixError::MalformedProfile);
        }
        scores.extend(
            fields
                .map(|field| {
                    field
                        .parse::<i32>()
                        .map_err(|_| MatrixError::MalformedProfile)
                })
                .collect::<Result<Vec<_>, _>>()?,
        );
    }
    if scores.len() == alphabet.len() * alphabet.len() {
        Ok((alphabet, scores))
    } else {
        Err(MatrixError::MalformedProfile)
    }
}

#[cfg(test)]
#[path = "profiles_tests.rs"]
mod tests;
