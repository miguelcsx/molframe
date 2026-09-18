//! Bounded network retrieval with mandatory content verification.

use sha2::{Digest as _, Sha256};
use std::io::Read as _;

/// Bytes retrieved only after their expected SHA-256 digest matched.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedDownload {
    /// Verified response body.
    pub bytes: Vec<u8>,
    /// Canonical lowercase digest.
    pub sha256: String,
}

/// Explicit resource and network policy for one verified download.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DownloadOptions {
    /// Maximum accepted decoded response body size.
    pub max_bytes: u64,
    /// Whole-request timeout.
    pub timeout: std::time::Duration,
    /// Maximum redirects; zero refuses redirects.
    pub redirect_limit: usize,
}

/// Failure while retrieving or verifying external bytes.
#[derive(Debug, thiserror::Error)]
pub enum DownloadError {
    /// URL is empty or the size bound is zero.
    #[error("download URL must be non-empty and max_bytes must be positive")]
    InvalidRequest,
    /// Expected checksum is not 64 hexadecimal digits.
    #[error("expected SHA-256 must contain exactly 64 hexadecimal digits")]
    InvalidChecksum,
    /// HTTP transport or status failed.
    #[error("HTTP download failed: {0}")]
    Http(#[from] reqwest::Error),
    /// Response body could not be read.
    #[error("download body failed: {0}")]
    Io(#[from] std::io::Error),
    /// Response exceeded the caller's bound.
    #[error("download exceeded {limit} bytes")]
    LimitExceeded {
        /// Explicit caller-provided byte bound.
        limit: u64,
    },
    /// Body did not match the pinned digest.
    #[error("SHA-256 mismatch: expected {expected}, received {actual}")]
    ChecksumMismatch {
        /// Pinned expected digest.
        expected: String,
        /// Digest of the received bytes.
        actual: String,
    },
}

/// Retrieves one URL with a byte ceiling and mandatory SHA-256 verification.
///
/// # Errors
///
/// Returns [`DownloadError`] for invalid controls, HTTP/I/O failure, size-limit
/// excess or checksum mismatch.
pub fn fetch_verified(
    url: &str,
    expected_sha256: &str,
    options: DownloadOptions,
) -> Result<VerifiedDownload, DownloadError> {
    if url.is_empty() || options.max_bytes == 0 || options.timeout.is_zero() {
        return Err(DownloadError::InvalidRequest);
    }
    let expected = normalise_checksum(expected_sha256)?;
    let redirect = if options.redirect_limit == 0 {
        reqwest::redirect::Policy::none()
    } else {
        reqwest::redirect::Policy::limited(options.redirect_limit)
    };
    let response = reqwest::blocking::Client::builder()
        .timeout(options.timeout)
        .redirect(redirect)
        .build()?
        .get(url)
        .send()?
        .error_for_status()?;
    if response
        .content_length()
        .is_some_and(|length| length > options.max_bytes)
    {
        return Err(DownloadError::LimitExceeded {
            limit: options.max_bytes,
        });
    }
    let read_limit = options
        .max_bytes
        .checked_add(1)
        .ok_or(DownloadError::InvalidRequest)?;
    let mut bytes = Vec::new();
    response.take(read_limit).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > options.max_bytes {
        return Err(DownloadError::LimitExceeded {
            limit: options.max_bytes,
        });
    }
    let actual = digest(&bytes);
    if actual != expected {
        return Err(DownloadError::ChecksumMismatch { expected, actual });
    }
    Ok(VerifiedDownload {
        bytes,
        sha256: actual,
    })
}

fn normalise_checksum(value: &str) -> Result<String, DownloadError> {
    let value = value.to_ascii_lowercase();
    if value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        Ok(value)
    } else {
        Err(DownloadError::InvalidChecksum)
    }
}

fn digest(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut output = String::with_capacity(64);
    for byte in digest {
        use std::fmt::Write as _;
        let _ = write!(output, "{byte:02x}");
    }
    output
}

#[cfg(test)]
#[path = "download_tests.rs"]
mod tests;
