//! Recognising a format, and the contract every reader and writer satisfies.
//!
//! Content decides before the file name does. A `.cif` suffix is used by the
//! archive format, by computed-model files, by small-molecule files and by
//! component definitions, so dispatching on it alone picks the wrong reader
//! regularly. A reader that can recognise its own content says so, and only when
//! nothing recognises the content does the suffix get a vote.
//!
//! Reading returns findings alongside success. Short-circuiting on the first
//! problem makes real archive files unusable; discarding findings on success
//! makes their problems invisible.

use super::input::{InputBuffer, Limits};
use crate::diagnostic::{Code, Diagnostic, Strictness};
use crate::structure::Structure;

/// A format pdbiox can read or write.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[non_exhaustive]
pub enum Format {
    /// Decide from the content, then from the file name.
    #[default]
    Auto,
    /// The archive's structured text format.
    Mmcif,
    /// The legacy fixed-column format.
    Pdb,
}

impl Format {
    /// The formats a caller may name, excluding automatic detection.
    pub const NAMED: [Self; 2] = [Self::Mmcif, Self::Pdb];

    /// The format's short name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Mmcif => "mmcif",
            Self::Pdb => "pdb",
        }
    }

    /// The file suffixes conventionally used for this format.
    #[must_use]
    pub const fn extensions(self) -> &'static [&'static str] {
        match self {
            Self::Auto => &[],
            Self::Mmcif => &["cif", "mmcif"],
            Self::Pdb => &["pdb", "ent"],
        }
    }

    /// Parses a format name.
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "auto" => Some(Self::Auto),
            "mmcif" | "cif" => Some(Self::Mmcif),
            "pdb" | "ent" => Some(Self::Pdb),
            _ => None,
        }
    }

    /// Returns true when this format's content is recognisable in `bytes`.
    ///
    /// Structured-text files open with a block header. Fixed-column files open
    /// with a record name in the first columns, and the coordinate records are
    /// the ones that must be there for the file to be worth reading.
    #[must_use]
    pub fn recognises(self, bytes: &[u8]) -> bool {
        match self {
            Self::Auto => false,
            Self::Mmcif => first_lines(bytes, 8).any(|line| line.starts_with(b"data_")),
            Self::Pdb => first_lines(bytes, 64).any(|line| {
                [&b"ATOM  "[..], b"HETATM", b"HEADER", b"CRYST1", b"MODEL "]
                    .iter()
                    .any(|record| line.starts_with(record))
            }),
        }
    }

    /// Decides the format of `input`, given what the caller asked for and the
    /// name the bytes came from.
    ///
    /// # Errors
    ///
    /// Returns a finding naming what was tried when nothing recognises the
    /// content and the name settles nothing.
    pub fn detect(
        requested: Self,
        input: &InputBuffer,
        name: Option<&str>,
    ) -> Result<Self, Diagnostic> {
        if requested != Self::Auto {
            return Ok(requested);
        }
        if let Some(format) = Self::NAMED
            .into_iter()
            .find(|f| f.recognises(input.as_bytes()))
        {
            return Ok(format);
        }
        if let Some(format) = name.and_then(Self::from_name) {
            return Ok(format);
        }
        Err(Diagnostic::new(Code::E1001)
            .with_context("tried", "content, then file name")
            .with_context(
                "name",
                match name {
                    Some(name) => name,
                    None => "(none)",
                },
            ))
    }

    /// The format a file name suggests, ignoring any compression suffix.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        let mut suffixes = name.rsplit('.');
        let mut candidate = suffixes.next()?;
        if matches!(candidate, "gz" | "zst" | "bz2" | "xz") {
            candidate = suffixes.next()?;
        }
        let lowered = candidate.to_ascii_lowercase();
        Self::NAMED
            .into_iter()
            .find(|format| format.extensions().contains(&lowered.as_str()))
    }
}

fn first_lines(bytes: &[u8], count: usize) -> impl Iterator<Item = &[u8]> {
    bytes.split(|byte| *byte == b'\n').take(count)
}

/// How much irregularity a read tolerates, and what it does about it.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
#[non_exhaustive]
pub enum ParseMode {
    /// Anything violating the specification fails the read.
    Strict,
    /// Produce a structure and report what was wrong with it.
    #[default]
    Permissive,
    /// Continue past local errors, marking the regions they affected.
    Recover,
}

impl ParseMode {
    /// The strictness this mode implies.
    #[must_use]
    pub const fn strictness(self) -> Strictness {
        match self {
            Self::Strict => Strictness::Strict,
            Self::Permissive => Strictness::Medium,
            Self::Recover => Strictness::Loose,
        }
    }
}

/// What a read should and should not bother doing.
///
/// Skipping work the caller does not want is the largest single lever on read
/// time. A caller who wants coordinates and nothing else should not pay for
/// citations, and says so here.
///
/// # Examples
///
/// ```
/// use pdbiox_core::io::{Format, ParseMode, ReadOptions};
///
/// let options = ReadOptions::new()
///     .format(Format::Mmcif)
///     .mode(ParseMode::Strict)
///     .only_first_model(true)
///     .only_atomic_coords(true);
///
/// assert_eq!(options.format, Format::Mmcif);
/// assert!(options.only_first_model);
/// ```
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    /// Which format, or automatic detection.
    pub format: Format,
    /// How much irregularity to tolerate.
    pub mode: ParseMode,
    /// Read only the first model.
    pub only_first_model: bool,
    /// Read coordinates and skip everything else.
    pub only_atomic_coords: bool,
    /// Drop hydrogens while reading.
    pub discard_hydrogens: bool,
    /// Ceilings that apply while reading.
    pub limits: Limits,
}

impl ReadOptions {
    /// Options that read everything, tolerating what the archive contains.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the format.
    #[must_use]
    pub const fn format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Sets how much irregularity to tolerate.
    #[must_use]
    pub const fn mode(mut self, mode: ParseMode) -> Self {
        self.mode = mode;
        self
    }

    /// Reads only the first model.
    #[must_use]
    pub const fn only_first_model(mut self, only: bool) -> Self {
        self.only_first_model = only;
        self
    }

    /// Reads coordinates and skips everything else.
    #[must_use]
    pub const fn only_atomic_coords(mut self, only: bool) -> Self {
        self.only_atomic_coords = only;
        self
    }

    /// Drops hydrogens while reading.
    #[must_use]
    pub const fn discard_hydrogens(mut self, discard: bool) -> Self {
        self.discard_hydrogens = discard;
        self
    }

    /// Sets the ceilings that apply while reading.
    #[must_use]
    pub const fn limits(mut self, limits: Limits) -> Self {
        self.limits = limits;
        self
    }
}

/// A structure and everything that was wrong with the file it came from, or the
/// findings that stopped it being read at all.
pub type ReadResult = Result<(Structure, Vec<Diagnostic>), Vec<Diagnostic>>;

/// A reader for one format.
pub trait Reader {
    /// The format this reads.
    const FORMAT: Format;

    /// Reads a structure.
    ///
    /// # Errors
    ///
    /// Returns the findings that stopped the read.
    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult;
}

/// Which parts of a structure a write should include.
///
/// Every method accepts by default, so a filter overrides only what it cares
/// about. Rejecting a parent excludes its children without the filter having to
/// say so twice.
pub trait Select {
    /// Whether to write this model.
    fn accept_model(&self, _model: crate::index::ModelIndex) -> bool {
        true
    }
    /// Whether to write this chain.
    fn accept_chain(&self, _chain: crate::index::ChainIndex) -> bool {
        true
    }
    /// Whether to write this residue.
    fn accept_residue(&self, _residue: crate::index::ResidueIndex) -> bool {
        true
    }
    /// Whether to write this atom.
    fn accept_atom(&self, _atom: crate::index::AtomIndex) -> bool {
        true
    }
}

/// A filter that writes everything.
#[derive(Clone, Copy, Debug, Default)]
pub struct SelectAll;

impl Select for SelectAll {}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
