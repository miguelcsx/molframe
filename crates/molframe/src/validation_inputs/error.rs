//! Why a validation reference input could not be loaded.

use crate::document::DocumentError;
use molframe_core::Findings;
use molframe_validate::{ReferenceError, RotamerError};

/// Failure while loading a validation reference input.
#[derive(Debug, thiserror::Error)]
pub enum ValidationInputError {
    /// The file could not be read.
    #[error("validation input I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// JSON syntax or schema was invalid.
    #[error("invalid JSON validation input: {0}")]
    Json(#[from] serde_json::Error),
    /// TOML syntax or schema was invalid.
    #[error("invalid TOML validation input: {0}")]
    Toml(#[from] toml::de::Error),
    /// The file suffix names neither JSON nor TOML.
    #[error("validation inputs must use a .json or .toml suffix")]
    UnsupportedFormat,
    /// A field held a value its schema does not allow.
    #[error("invalid value for {field}: {reason}")]
    InvalidValue {
        /// Dotted field path within the document, such as `distribution.edges`.
        field: String,
        /// What is wrong with the value.
        reason: String,
    },
    /// A reference library or one of its distributions was refused.
    #[error("invalid reference library: {0}")]
    Reference(#[from] ReferenceError),
    /// A rotamer profile was refused.
    #[error("invalid rotamer profile: {0}")]
    Rotamer(#[from] RotamerError),
    /// A selection expression failed to compile or evaluate.
    #[error("selection for {id} failed: {findings}")]
    Selection {
        /// The plane or group the selection belongs to.
        id: String,
        /// The compile or evaluation diagnostics.
        findings: Findings,
    },
}

impl From<DocumentError> for ValidationInputError {
    fn from(error: DocumentError) -> Self {
        match error {
            DocumentError::Io(error) => Self::Io(error),
            DocumentError::Json(error) => Self::Json(error),
            DocumentError::Toml(error) => Self::Toml(error),
            DocumentError::UnsupportedFormat => Self::UnsupportedFormat,
        }
    }
}

impl ValidationInputError {
    pub(super) fn invalid(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidValue {
            field: field.into(),
            reason: reason.into(),
        }
    }
}
