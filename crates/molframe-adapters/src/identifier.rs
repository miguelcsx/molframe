//! Identifier acquisition boundary.
//!
//! A direct URL can be retrieved by [`crate::fetch_verified`], but this
//! workspace does not contain an authoritative identifier-to-resource index
//! with pinned digests. Keeping that distinction explicit prevents a caller
//! from mistaking a guessed URL or a locally computed digest for verified
//! acquisition.

use crate::VerifiedDownload;

/// Failure while resolving a public structure identifier.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum IdentifierFetchError {
    /// The caller supplied no identifier.
    #[error("structure identifier must be non-empty")]
    EmptyIdentifier,
    /// No authoritative provider manifest is configured for this build.
    #[error("no authoritative provider is configured for structure identifier {identifier:?}")]
    ProviderUnavailable {
        /// Identifier that could not be resolved.
        identifier: String,
    },
}

/// Resolve and fetch a structure by identifier.
///
/// This deliberately reports [`IdentifierFetchError::ProviderUnavailable`]
/// until a provider supplies both a canonical URL and an authoritative
/// checksum/resource index. It does not guess a URL, infer a format, or treat
/// a digest computed after download as publisher authentication. Callers with
/// an already verified URL should use [`crate::fetch_verified`] instead.
///
/// # Errors
///
/// Returns [`IdentifierFetchError::EmptyIdentifier`] for an empty identifier
/// and [`IdentifierFetchError::ProviderUnavailable`] for every non-empty
/// identifier while no provider manifest is configured.
pub fn fetch_identifier(identifier: &str) -> Result<VerifiedDownload, IdentifierFetchError> {
    if identifier.trim().is_empty() {
        return Err(IdentifierFetchError::EmptyIdentifier);
    }
    Err(IdentifierFetchError::ProviderUnavailable {
        identifier: identifier.to_owned(),
    })
}

#[cfg(test)]
#[path = "identifier_tests.rs"]
mod tests;
