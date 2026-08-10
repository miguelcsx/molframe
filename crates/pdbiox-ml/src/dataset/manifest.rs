//! Manifest parsing and lazy entry access.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use super::{DatasetFilter, DatasetSplit, SplitOptions};

/// One structure handle and its load-free statistics.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ManifestEntry {
    /// Stable unique identifier.
    pub id: Box<str>,
    /// Structure path or resolver key.
    pub path: PathBuf,
    /// Atom count available without loading the structure.
    pub atom_count: u64,
    /// Experimental resolution in ångström, when applicable.
    #[serde(default)]
    pub resolution: Option<f32>,
    /// Experimental or prediction method.
    #[serde(default)]
    pub method: Option<Box<str>>,
    /// ISO-8601 deposition date (`YYYY-MM-DD`).
    #[serde(default)]
    pub deposition_date: Option<Box<str>>,
    /// Sequence used by identity-aware splitting.
    #[serde(default)]
    pub sequence: Option<Box<str>>,
    /// Precomputed structural cluster identifier.
    #[serde(default)]
    pub structure_cluster: Option<Box<str>>,
    /// User-defined tags.
    #[serde(default)]
    pub tags: Vec<Box<str>>,
    /// Numeric statistics used by downstream manifest queries.
    #[serde(default)]
    pub statistics: BTreeMap<Box<str>, f64>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum ManifestDocument {
    Entries { entries: Vec<ManifestEntry> },
    Array(Vec<ManifestEntry>),
}

/// Immutable lazy dataset backed only by manifest entries and selected indices.
#[derive(Clone, Debug, PartialEq)]
pub struct Dataset {
    entries: Arc<[ManifestEntry]>,
    indices: Arc<[usize]>,
}

impl Dataset {
    /// Builds a dataset after validating identifiers, paths and statistics.
    ///
    /// # Errors
    ///
    /// Rejects duplicate or empty identifiers, empty paths, invalid resolution,
    /// malformed dates and non-finite statistics.
    pub fn new(entries: Vec<ManifestEntry>) -> Result<Self, DatasetError> {
        validate_entries(&entries)?;
        let indices = (0..entries.len()).collect::<Vec<_>>();
        Ok(Self {
            entries: entries.into(),
            indices: indices.into(),
        })
    }

    /// Reads a JSON manifest without loading any structure.
    ///
    /// # Errors
    ///
    /// Returns filesystem, JSON-schema or manifest validation errors.
    pub fn from_manifest(path: &Path) -> Result<Self, DatasetError> {
        let bytes = std::fs::read(path)?;
        let document: ManifestDocument = serde_json::from_slice(&bytes)?;
        let mut entries = match document {
            ManifestDocument::Entries { entries } | ManifestDocument::Array(entries) => entries,
        };
        if let Some(parent) = path.parent() {
            for entry in &mut entries {
                if entry.path.is_relative() {
                    entry.path = parent.join(&entry.path);
                }
            }
        }
        Self::new(entries)
    }

    /// Number of selected manifest entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Whether no entries are selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Selected entries in deterministic manifest order.
    #[must_use]
    pub fn entries(&self) -> impl ExactSizeIterator<Item = &ManifestEntry> {
        self.indices.iter().map(|index| &self.entries[*index])
    }

    /// Returns a filtered view without reading coordinates.
    ///
    /// # Errors
    ///
    /// Rejects non-positive/non-finite resolution bounds and inverted atom
    /// count ranges.
    pub fn filter(&self, filter: &DatasetFilter) -> Result<Self, DatasetError> {
        filter.validate()?;
        Ok(self.subset(
            self.indices
                .iter()
                .copied()
                .filter(|index| filter.matches(&self.entries[*index]))
                .collect(),
        ))
    }

    /// Loads one selected entry through the caller's format resolver.
    ///
    /// # Errors
    ///
    /// Propagates the loader error or returns [`DatasetError::IndexOutOfBounds`].
    pub fn load_with<T, E>(
        &self,
        index: usize,
        loader: impl FnOnce(&ManifestEntry) -> Result<T, E>,
    ) -> Result<T, LoadError<E>> {
        let absolute = self
            .indices
            .get(index)
            .copied()
            .ok_or(LoadError::Dataset(DatasetError::IndexOutOfBounds { index }))?;
        loader(&self.entries[absolute]).map_err(LoadError::Loader)
    }

    /// Splits the manifest without loading coordinates.
    ///
    /// # Errors
    ///
    /// Rejects invalid ratios or metadata missing for the selected strategy.
    pub fn split(&self, options: &SplitOptions) -> Result<DatasetSplit, DatasetError> {
        super::split::split(self, options)
    }

    /// Divides selected handles into lazy batches without loading coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a zero batch size.
    pub fn batches(&self, batch_size: usize) -> Result<Vec<Self>, DatasetError> {
        if batch_size == 0 {
            return Err(DatasetError::InvalidBatchSize);
        }
        Ok(self
            .indices
            .chunks(batch_size)
            .map(|indices| self.subset(indices.to_vec()))
            .collect())
    }

    pub(crate) fn absolute_indices(&self) -> &[usize] {
        &self.indices
    }

    pub(crate) fn entry(&self, index: usize) -> &ManifestEntry {
        &self.entries[index]
    }

    pub(crate) fn subset(&self, indices: Vec<usize>) -> Self {
        Self {
            entries: self.entries.clone(),
            indices: indices.into(),
        }
    }
}

/// Dataset or caller loader failure.
#[derive(Debug, thiserror::Error)]
pub enum LoadError<E> {
    /// Dataset index failure.
    #[error(transparent)]
    Dataset(#[from] DatasetError),
    /// Caller-provided loader failure.
    #[error("dataset loader failed: {0}")]
    Loader(E),
}

/// Manifest and splitting errors.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum DatasetError {
    /// Exact sequence alignment exceeded its supported numeric domain.
    #[error(transparent)]
    Alignment(#[from] pdbiox_seq::AlignError),
    /// Manifest filesystem failure.
    #[error("dataset manifest I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// Manifest JSON failure.
    #[error("invalid dataset manifest JSON: {0}")]
    Json(#[from] serde_json::Error),
    /// An entry violates the manifest contract.
    #[error("invalid dataset entry at index {index}: {reason}")]
    InvalidEntry {
        /// Manifest index.
        index: usize,
        /// Stable reason.
        reason: &'static str,
    },
    /// Stable identifiers must be unique.
    #[error("duplicate dataset identifier: {id}")]
    DuplicateId {
        /// Repeated identifier.
        id: Box<str>,
    },
    /// Selected index is absent.
    #[error("dataset index {index} is out of bounds")]
    IndexOutOfBounds {
        /// Requested selected index.
        index: usize,
    },
    /// Split ratios do not form a valid distribution.
    #[error("invalid dataset split ratios")]
    InvalidRatios,
    /// A split strategy requires absent metadata.
    #[error("entry {id} lacks {field} required by the split strategy")]
    MissingSplitMetadata {
        /// Entry identifier.
        id: Box<str>,
        /// Required field.
        field: &'static str,
    },
    /// Sequence threshold is outside the closed unit interval.
    #[error("sequence identity threshold must be between zero and one")]
    InvalidThreshold,
    /// A lazy batch must contain at least one handle.
    #[error("dataset batch size must be greater than zero")]
    InvalidBatchSize,
    /// A filter contains invalid numeric bounds.
    #[error("invalid dataset filter")]
    InvalidFilter,
}

fn validate_entries(entries: &[ManifestEntry]) -> Result<(), DatasetError> {
    let mut identifiers = std::collections::BTreeSet::new();
    for (index, entry) in entries.iter().enumerate() {
        let invalid = if entry.id.trim().is_empty() {
            Some("identifier is empty")
        } else if entry.path.as_os_str().is_empty() {
            Some("path is empty")
        } else if entry.atom_count == 0 {
            Some("atom count is zero")
        } else if entry
            .resolution
            .is_some_and(|value| !value.is_finite() || value <= 0.0)
        {
            Some("resolution is not positive and finite")
        } else if entry
            .sequence
            .as_deref()
            .is_some_and(|sequence| sequence.is_empty() || !sequence.is_ascii())
        {
            Some("sequence is empty or non-ASCII")
        } else if entry
            .deposition_date
            .as_deref()
            .is_some_and(|date| !valid_date(date))
        {
            Some("deposition date is not YYYY-MM-DD")
        } else if entry.statistics.values().any(|value| !value.is_finite()) {
            Some("statistic is not finite")
        } else {
            None
        };
        if let Some(reason) = invalid {
            return Err(DatasetError::InvalidEntry { index, reason });
        }
        if !identifiers.insert(entry.id.as_ref()) {
            return Err(DatasetError::DuplicateId {
                id: entry.id.clone(),
            });
        }
    }
    Ok(())
}

fn valid_date(date: &str) -> bool {
    let bytes = date.as_bytes();
    if bytes.len() != 10 || bytes[4] != b'-' || bytes[7] != b'-' {
        return false;
    }
    let numeric = bytes
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != 4 && *index != 7)
        .all(|(_, value)| value.is_ascii_digit());
    if !numeric {
        return false;
    }
    let year = date[..4].parse::<u32>().ok();
    let month = date[5..7].parse::<u8>().ok();
    let day = date[8..10].parse::<u8>().ok();
    year.zip(month)
        .zip(day)
        .is_some_and(|((year, month), day)| day > 0 && day <= days_in_month(year, month))
}

const fn days_in_month(year: u32, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if year.is_multiple_of(400) || (year.is_multiple_of(4) && !year.is_multiple_of(100)) => {
            29
        }
        2 => 28,
        _ => 0,
    }
}
