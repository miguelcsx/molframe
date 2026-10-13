//! Converters from a downloaded file to the rows of a table.
//!
//! A converter is a named, deterministic text transformation: the same bytes
//! always give the same rows. `tsv` keeps tab-separated rows as they are. The
//! other ids live in their own modules. Cost: one pass over the text.

use super::stereo_chemical_props;

/// A converter id from the manifest.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Converter {
    /// Rows that are already in the table's layout.
    Tsv,
    /// A stereochemical parameter file of bond and angle sections.
    StereoChemicalProps,
}

impl Converter {
    /// The converter called `id`.
    pub(super) fn named(id: &str) -> Result<Self, String> {
        match id {
            "tsv" => Ok(Self::Tsv),
            "stereo_chemical_props" => Ok(Self::StereoChemicalProps),
            other => Err(format!(
                "unknown converter `{other}`; use tsv or stereo_chemical_props"
            )),
        }
    }

    /// The rows `text` converts to, each without its line ending.
    pub(super) fn rows(self, text: &str) -> Result<Vec<String>, String> {
        match self {
            Self::Tsv => Ok(tab_rows(text)),
            Self::StereoChemicalProps => stereo_chemical_props::rows(text),
        }
    }
}

fn tab_rows(text: &str) -> Vec<String> {
    text.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

#[cfg(test)]
#[path = "convert_tests.rs"]
mod tests;
