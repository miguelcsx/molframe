//! Policies for reusable scratch memory and deterministic spill storage.

use std::path::{Path, PathBuf};

/// Upper bound for reusable in-memory scratch owned by one stage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScratchPolicy {
    max_bytes: usize,
}

impl ScratchPolicy {
    /// Creates a scratch policy bounded by `max_bytes`.
    #[must_use]
    pub const fn new(max_bytes: usize) -> Self {
        Self { max_bytes }
    }

    /// Returns the maximum reusable scratch capacity.
    #[must_use]
    pub const fn max_bytes(self) -> usize {
        self.max_bytes
    }
}

impl Default for ScratchPolicy {
    fn default() -> Self {
        Self::new(8_000_000)
    }
}

/// Policy controlling whether and where execution may spill temporary data.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum TempStoragePolicy {
    /// Operations that cannot remain within memory must return an error.
    #[default]
    Disabled,
    /// Spill files may be created below one explicit directory.
    Directory {
        /// Root below which implementations create execution-specific files.
        root: PathBuf,
        /// Maximum bytes retained on disk by one execution context.
        max_bytes: u64,
    },
}

impl TempStoragePolicy {
    /// Enables bounded spill below `root`.
    #[must_use]
    pub fn directory(root: impl Into<PathBuf>, max_bytes: u64) -> Self {
        Self::Directory {
            root: root.into(),
            max_bytes,
        }
    }

    /// Returns the configured spill root, when disk storage is enabled.
    #[must_use]
    pub fn root(&self) -> Option<&Path> {
        match self {
            Self::Disabled => None,
            Self::Directory { root, .. } => Some(root),
        }
    }

    /// Returns the disk-space ceiling. Disabled storage has a zero ceiling.
    #[must_use]
    pub const fn max_bytes(&self) -> u64 {
        match self {
            Self::Disabled => 0,
            Self::Directory { max_bytes, .. } => *max_bytes,
        }
    }
}

#[cfg(test)]
#[path = "storage_tests.rs"]
mod tests;
