//! Content and suffix based format recognition.

use super::super::input::InputBuffer;
use crate::diagnostic::{Code, Diagnostic};

const PDBML_PROBE_WINDOWS: usize = 2_048;
const MESSAGEPACK_PROBE_WINDOWS: usize = 512;

/// A format pdbiox can read or write.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[non_exhaustive]
pub enum Format {
    /// Decide from the content, then from the file name.
    #[default]
    Auto,
    /// The archive's structured text format.
    Mmcif,
    /// XML serialization of the `PDBx` data model.
    Pdbml,
    /// The archive's `MessagePack`-based binary format.
    BinaryCif,
    /// The deprecated but still encountered macromolecular transmission format.
    Mmtf,
    /// The legacy fixed-column format.
    Pdb,
    /// PDB-shaped coordinates carrying partial charge and radius.
    Pqr,
    /// `AutoDock`'s PDB-shaped coordinates carrying charge and atom type.
    Pdbqt,
}

impl Format {
    /// The formats a caller may name, excluding automatic detection.
    pub const NAMED: [Self; 7] = [
        Self::Mmcif,
        Self::Pdbml,
        Self::BinaryCif,
        Self::Mmtf,
        Self::Pdbqt,
        Self::Pqr,
        Self::Pdb,
    ];

    /// The format's short name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Mmcif => "mmcif",
            Self::Pdbml => "pdbml",
            Self::BinaryCif => "bcif",
            Self::Mmtf => "mmtf",
            Self::Pdb => "pdb",
            Self::Pqr => "pqr",
            Self::Pdbqt => "pdbqt",
        }
    }

    /// The file suffixes conventionally used for this format.
    #[must_use]
    pub const fn extensions(self) -> &'static [&'static str] {
        match self {
            Self::Auto => &[],
            Self::Mmcif => &["cif", "mmcif"],
            Self::Pdbml => &["xml", "pdbml"],
            Self::BinaryCif => &["bcif"],
            Self::Mmtf => &["mmtf"],
            Self::Pdb => &["pdb", "ent"],
            Self::Pqr => &["pqr"],
            Self::Pdbqt => &["pdbqt"],
        }
    }

    /// Parses a format name without allocating a normalized copy.
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        if name.eq_ignore_ascii_case("auto") {
            Some(Self::Auto)
        } else if name.eq_ignore_ascii_case("mmcif") || name.eq_ignore_ascii_case("cif") {
            Some(Self::Mmcif)
        } else if name.eq_ignore_ascii_case("pdbml") || name.eq_ignore_ascii_case("xml") {
            Some(Self::Pdbml)
        } else if name.eq_ignore_ascii_case("bcif") || name.eq_ignore_ascii_case("binarycif") {
            Some(Self::BinaryCif)
        } else if name.eq_ignore_ascii_case("mmtf") {
            Some(Self::Mmtf)
        } else if name.eq_ignore_ascii_case("pdb") || name.eq_ignore_ascii_case("ent") {
            Some(Self::Pdb)
        } else if name.eq_ignore_ascii_case("pqr") {
            Some(Self::Pqr)
        } else if name.eq_ignore_ascii_case("pdbqt") {
            Some(Self::Pdbqt)
        } else {
            None
        }
    }

    /// Returns true when this format's content is recognisable in `bytes`.
    #[must_use]
    pub fn recognises(self, bytes: &[u8]) -> bool {
        match self {
            Self::Auto | Self::Pqr => false,
            Self::Mmcif => first_lines(bytes, 8).any(|line| line.starts_with(b"data_")),
            Self::Pdbml => contains_in_probe(bytes, b"PDBx:datablock", PDBML_PROBE_WINDOWS),
            Self::BinaryCif => recognises_messagepack_map(bytes, b"dataBlocks"),
            Self::Mmtf => recognises_messagepack_map(bytes, b"mmtfVersion"),
            Self::Pdbqt => recognises_pdbqt(bytes),
            Self::Pdb => recognises_pdb(bytes),
        }
    }

    /// Decides the format of `input`, given what the caller asked for and the
    /// name the bytes came from.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when neither content nor the optional name identifies
    /// a supported format.
    pub fn detect(
        requested: Self,
        input: &InputBuffer,
        name: Option<&str>,
    ) -> Result<Self, Diagnostic> {
        if requested != Self::Auto {
            return Ok(requested);
        }

        let bytes = input.as_bytes();
        for format in [
            Self::Mmcif,
            Self::Pdbml,
            Self::BinaryCif,
            Self::Mmtf,
            Self::Pdbqt,
        ] {
            if format.recognises(bytes) {
                return Ok(format);
            }
        }

        let named = name.and_then(Self::from_name);
        if recognises_pdb(bytes) {
            if let Some(format @ (Self::Pqr | Self::Pdbqt)) = named {
                return Ok(format);
            }
            return Ok(Self::Pdb);
        }
        if let Some(format) = named {
            return Ok(format);
        }

        let name = match name {
            Some(name) => name,
            None => "(none)",
        };
        Err(Diagnostic::new(Code::E1001)
            .with_context("tried", "content, then file name")
            .with_context("name", name))
    }

    /// The format a file name suggests, ignoring any compression suffix.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        let mut suffixes = name.rsplit('.');
        let mut candidate = suffixes.next()?;
        if is_compression_suffix(candidate) {
            candidate = suffixes.next()?;
        }
        Self::NAMED.into_iter().find(|format| {
            format
                .extensions()
                .iter()
                .any(|extension| candidate.eq_ignore_ascii_case(extension))
        })
    }
}

fn is_compression_suffix(suffix: &str) -> bool {
    suffix.eq_ignore_ascii_case("gz")
        || suffix.eq_ignore_ascii_case("zst")
        || suffix.eq_ignore_ascii_case("bz2")
        || suffix.eq_ignore_ascii_case("xz")
}

fn is_messagepack_map_prefix(byte: u8) -> bool {
    (0x80..=0x8f).contains(&byte) || matches!(byte, 0xde | 0xdf)
}

fn recognises_messagepack_map(bytes: &[u8], key: &[u8]) -> bool {
    bytes
        .first()
        .copied()
        .is_some_and(is_messagepack_map_prefix)
        && contains_in_probe(bytes, key, MESSAGEPACK_PROBE_WINDOWS)
}

fn contains_in_probe(bytes: &[u8], needle: &[u8], windows: usize) -> bool {
    if needle.is_empty() {
        return true;
    }
    bytes
        .windows(needle.len())
        .take(windows)
        .any(|window| window == needle)
}

fn recognises_pdb(bytes: &[u8]) -> bool {
    first_lines(bytes, 64).any(|line| {
        line.starts_with(b"ATOM  ")
            || line.starts_with(b"HETATM")
            || line.starts_with(b"HEADER")
            || line.starts_with(b"CRYST1")
            || line.starts_with(b"MODEL ")
    })
}

fn recognises_pdbqt(bytes: &[u8]) -> bool {
    first_lines(bytes, 64).any(recognises_pdbqt_line)
}

fn recognises_pdbqt_line(line: &[u8]) -> bool {
    if matches!(line, b"ROOT" | b"ENDROOT") || line.starts_with(b"TORSDOF") {
        return true;
    }
    if !(line.starts_with(b"ATOM  ") || line.starts_with(b"HETATM")) {
        return false;
    }
    let charge = line
        .get(70..76)
        .and_then(|field| str::from_utf8(field).ok());
    let atom_type = line.get(77..).and_then(|field| str::from_utf8(field).ok());
    charge.is_some_and(|field| field.trim().parse::<f64>().is_ok())
        && atom_type.is_some_and(|field| !field.trim().is_empty())
}

fn first_lines(bytes: &[u8], count: usize) -> impl Iterator<Item = &[u8]> {
    bytes
        .split(|byte| *byte == b'\n')
        .take(count)
        .map(strip_carriage_return)
}

fn strip_carriage_return(line: &[u8]) -> &[u8] {
    if line.last().copied() == Some(b'\r') {
        &line[..line.len() - 1]
    } else {
        line
    }
}
