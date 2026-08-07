//! Bounded amortized-linear byte collection and expansion checks.

use super::Limits;
use crate::diagnostic::{Code, Diagnostic};
use std::collections::TryReserveError;
use std::io::{ErrorKind, Read};
use std::path::Path;
use std::sync::Arc;

const READ_BUFFER_BYTES: usize = 64 * 1024;

#[derive(Clone, Copy)]
pub(super) enum OutputLimit {
    Uncompressed {
        maximum: u64,
    },
    #[cfg(any(feature = "gzip", feature = "zstd"))]
    Expanded {
        maximum: u64,
        compressed: u64,
        limits: Limits,
    },
}

impl OutputLimit {
    pub(super) const fn uncompressed(limits: Limits) -> Self {
        Self::Uncompressed {
            maximum: limits.decompressed_bytes,
        }
    }

    #[cfg(any(feature = "gzip", feature = "zstd"))]
    pub(super) fn expanded(compressed: u64, limits: Limits) -> Self {
        Self::Expanded {
            maximum: maximum_expanded_bytes(compressed, limits),
            compressed,
            limits,
        }
    }

    const fn maximum(self) -> u64 {
        match self {
            Self::Uncompressed { maximum } => maximum,
            #[cfg(any(feature = "gzip", feature = "zstd"))]
            Self::Expanded { maximum, .. } => maximum,
        }
    }

    fn validate(self, expanded: u64) -> Result<(), Diagnostic> {
        match self {
            Self::Uncompressed { maximum } if expanded > maximum => {
                Err(Limits::exceeded("decompressed bytes", expanded))
            }
            Self::Uncompressed { .. } => Ok(()),
            #[cfg(any(feature = "gzip", feature = "zstd"))]
            Self::Expanded {
                compressed, limits, ..
            } => check_expansion(expanded, compressed, limits),
        }
    }
}

pub(super) fn collect_checked<R: Read>(
    mut reader: R,
    limit: OutputLimit,
    failure_message: &'static str,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    let mut output = Vec::new();
    let mut buffer = vec![0_u8; READ_BUFFER_BYTES].into_boxed_slice();
    loop {
        let Ok(current) = u64::try_from(output.len()) else {
            return Err(Limits::exceeded("decompressed bytes", output.len()));
        };
        let remaining = limit.maximum().saturating_sub(current);
        let request = match usize::try_from(remaining) {
            Ok(value) if value < buffer.len() => value + 1,
            Ok(_) | Err(_) => buffer.len(),
        };
        let count = match reader.read(&mut buffer[..request]) {
            Ok(0) => break,
            Ok(count) => count,
            Err(error) if error.kind() == ErrorKind::Interrupted => continue,
            Err(error) => return Err(read_failure(failure_message, path, &error)),
        };
        let Ok(additional) = u64::try_from(count) else {
            return Err(Limits::exceeded("decompressed bytes", count));
        };
        let Some(expanded) = current.checked_add(additional) else {
            return Err(Limits::exceeded("decompressed bytes", "more than u64::MAX"));
        };
        limit.validate(expanded)?;
        reserve_output(&mut output, count, limit.maximum(), path)?;
        output.extend_from_slice(&buffer[..count]);
    }
    Ok(output)
}

fn reserve_output(
    output: &mut Vec<u8>,
    additional: usize,
    maximum: u64,
    path: Option<&Path>,
) -> Result<(), Diagnostic> {
    let Some(required) = output.len().checked_add(additional) else {
        return Err(Limits::exceeded(
            "decompressed bytes",
            "more than usize::MAX",
        ));
    };
    if required <= output.capacity() {
        return Ok(());
    }
    let maximum_capacity = match usize::try_from(maximum) {
        Ok(value) => value,
        Err(_) => usize::MAX,
    };
    if required > maximum_capacity {
        return Err(Limits::exceeded("decompressed bytes", required));
    }
    let grown = output.capacity().saturating_mul(2).max(required);
    let target = grown.min(maximum_capacity);
    let Some(extra) = target.checked_sub(output.len()) else {
        return Err(Limits::exceeded("decompressed bytes", required));
    };
    output
        .try_reserve_exact(extra)
        .map_err(|error| allocation_failure(path, &error))
}

pub(super) fn read_failure(
    message: &'static str,
    path: Option<&Path>,
    error: &std::io::Error,
) -> Diagnostic {
    let finding = Diagnostic::new(Code::E1901).with_message(message);
    let finding = match path {
        Some(path) => finding.with_context("path", path.display().to_string()),
        None => finding,
    };
    finding.with_context("reason", error.to_string())
}

fn allocation_failure(path: Option<&Path>, error: &TryReserveError) -> Diagnostic {
    let finding = Diagnostic::new(Code::E1901).with_message("input buffer could not be allocated");
    let finding = match path {
        Some(path) => finding.with_context("path", path.display().to_string()),
        None => finding,
    };
    finding.with_context("reason", error.to_string())
}

pub(super) fn path_origin(path: &Path) -> Arc<str> {
    match path.to_str() {
        Some(origin) => Arc::from(origin),
        None => Arc::from(path.to_string_lossy().as_ref()),
    }
}

#[cfg(any(feature = "gzip", feature = "zstd"))]
fn maximum_expanded_bytes(compressed: u64, limits: Limits) -> u64 {
    limits
        .decompressed_bytes
        .min(ratio_byte_limit(compressed, limits.compression_ratio))
}

#[cfg(any(feature = "gzip", feature = "zstd", test))]
fn ratio_byte_limit(compressed: u64, ratio: u64) -> u64 {
    if compressed == 0 {
        return u64::MAX;
    }
    match compressed.checked_mul(ratio) {
        Some(value) => value,
        None => u64::MAX,
    }
}

#[cfg(any(feature = "gzip", feature = "zstd", test))]
pub(super) fn check_expansion(
    expanded: u64,
    compressed: u64,
    limits: Limits,
) -> Result<(), Diagnostic> {
    if expanded > limits.decompressed_bytes {
        return Err(Limits::exceeded("decompressed bytes", expanded));
    }
    if expanded > ratio_byte_limit(compressed, limits.compression_ratio) {
        return Err(Limits::exceeded(
            "compression ratio",
            ceiling_ratio(expanded, compressed),
        ));
    }
    Ok(())
}

#[cfg(any(feature = "gzip", feature = "zstd", test))]
fn ceiling_ratio(numerator: u64, denominator: u64) -> u64 {
    if denominator == 0 {
        return u64::MAX;
    }
    let quotient = numerator / denominator;
    if numerator.is_multiple_of(denominator) {
        quotient
    } else {
        quotient.saturating_add(1)
    }
}
