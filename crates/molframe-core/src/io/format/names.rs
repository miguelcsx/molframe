//! The words a read is configured with.

use super::detect::Format;
use super::options::{AmbiguousResidueBoundaryPolicy, MissingElementPolicy, ParseMode};
use crate::contract::{PolicyParseError, closed_vocabulary};
use std::fmt;
use std::str::FromStr;

closed_vocabulary!(
    ParseMode,
    "mode",
    [
        (ParseMode::Strict, "strict"),
        (ParseMode::Permissive, "permissive"),
        (ParseMode::Recover, "recover"),
    ]
);

closed_vocabulary!(
    MissingElementPolicy,
    "missing_element",
    [
        (MissingElementPolicy::PreserveUnknown, "preserve-unknown"),
        (
            MissingElementPolicy::InferFromAtomName,
            "infer-from-atom-name"
        ),
    ]
);

closed_vocabulary!(
    AmbiguousResidueBoundaryPolicy,
    "ambiguous_residue_boundary",
    [
        (AmbiguousResidueBoundaryPolicy::Reject, "reject"),
        (
            AmbiguousResidueBoundaryPolicy::InferFromFileOrder,
            "infer-from-file-order"
        ),
    ]
);

impl Format {
    /// The names a format can be spelled with, besides the file suffixes
    /// (`cif`, `ent`, `mol`, ...) that [`Format::parse`] also accepts.
    pub const NAMES: &'static [&'static str] = &[
        "auto", "mmcif", "pdbml", "bcif", "mmtf", "pdb", "pqr", "pdbqt", "sdf", "mol2", "smallcif",
    ];
}

impl fmt::Display for Format {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl FromStr for Format {
    type Err = PolicyParseError;

    /// Case-insensitive; a conventional file suffix names its format.
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::parse(value)
            .ok_or_else(|| PolicyParseError::new("format", value, &Self::NAMES.join(", ")))
    }
}

#[cfg(test)]
#[path = "names_tests.rs"]
mod tests;
