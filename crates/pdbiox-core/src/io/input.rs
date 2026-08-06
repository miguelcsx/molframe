//! Where bytes come from, and the limits that keep a hostile file from
//! exhausting the machine.
//!
//! Memory mapping is right for a large uncompressed file on local disk and wrong
//! for a compressed one, a network stream or bytes already in hand. All four are
//! ordinary here and none is architecturally privileged: a reader takes bytes and
//! does not care how they were obtained.
//!
//! Decompression is transparent and bounded. A file that expands a thousandfold
//! is refused rather than obeyed.

use crate::diagnostic::{Code, Diagnostic};
use std::collections::TryReserveError;
use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::Path;
use std::sync::Arc;

/// Default maximum accepted decompressed size: four gibibytes.
const DEFAULT_DECOMPRESSED_BYTES: u64 = 4_u64 << 30;

/// Default maximum decompressed-to-compressed size ratio.
const DEFAULT_COMPRESSION_RATIO: u64 = 1_000;

/// Default maximum declared row count for one category.
const DEFAULT_ROWS_PER_CATEGORY: u64 = 100_000_000;

/// Default maximum parser nesting depth.
const DEFAULT_NESTING_DEPTH: u32 = 64;

/// Default maximum identifier dictionary size.
const DEFAULT_DICTIONARY_ENTRIES: u32 = 1_000_000;

/// Stack buffer used while collecting input or decoder output.
const READ_BUFFER_BYTES: usize = 64 * 1024;

/// Number of leading bytes needed to identify every supported container.
const COMPRESSION_PREFIX_BYTES: usize = 4;

/// Ceilings that apply while reading.
///
/// These exist to turn a hostile input into a diagnostic rather than an
/// out-of-memory kill. The defaults are far above any real file.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Largest accepted size after decompression, in bytes.
    pub decompressed_bytes: u64,
    /// Largest accepted ratio of decompressed to compressed size.
    pub compression_ratio: u64,
    /// Largest accepted declared row count for one category.
    pub rows_per_category: u64,
    /// Deepest accepted nesting.
    pub nesting_depth: u32,
    /// Largest accepted identifier dictionary.
    pub dictionary_entries: u32,
}

impl Default for Limits {
    /// Returns conservative input limits suitable for real structure files.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    fn default() -> Self {
        Self {
            decompressed_bytes: DEFAULT_DECOMPRESSED_BYTES,
            compression_ratio: DEFAULT_COMPRESSION_RATIO,
            rows_per_category: DEFAULT_ROWS_PER_CATEGORY,
            nesting_depth: DEFAULT_NESTING_DEPTH,
            dictionary_entries: DEFAULT_DICTIONARY_ENTRIES,
        }
    }
}

impl Limits {
    /// The finding raised when a limit is exceeded.
    ///
    /// `limit` names the violated ceiling and `value` records the observed
    /// value. Formatting takes `O(value)` time and allocates only diagnostic
    /// context storage.
    #[must_use]
    pub fn exceeded(limit: &'static str, value: impl std::fmt::Display) -> Diagnostic {
        Diagnostic::new(Code::E1901)
            .with_context("limit", limit)
            .with_context("value", value.to_string())
    }
}

/// How a stream of bytes was compressed, where it was.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Compression {
    /// Not compressed.
    None,
    /// Gzip, as in a `.gz` suffix.
    Gzip,
    /// Zstandard, as in a `.zst` suffix.
    Zstd,
}

impl Compression {
    /// Recognises a compression container from its leading bytes.
    ///
    /// Content decides, not the file name: an archive mirror serving gzip under
    /// a plain `.cif` name is common enough to matter.
    ///
    /// Reads at most four bytes, runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn sniff(bytes: &[u8]) -> Self {
        match bytes {
            [0x1f, 0x8b, ..] => Self::Gzip,
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Self::Zstd,
            _ => Self::None,
        }
    }
}

/// Bytes to read, however they were obtained.
///
/// # Examples
///
/// ```
/// use pdbiox_core::io::InputBuffer;
///
/// let input = InputBuffer::from_bytes(b"data_test\n".to_vec());
/// assert_eq!(input.len(), 10);
/// assert!(input.as_bytes().starts_with(b"data_"));
/// ```
#[derive(Clone, Debug)]
pub struct InputBuffer {
    bytes: Arc<Vec<u8>>,
    origin: Option<Arc<str>>,
}

impl InputBuffer {
    /// Takes ownership of bytes already in hand.
    ///
    /// The original `Vec` allocation is retained without copying its contents.
    /// Runs in `O(1)` time and allocates only the reference-counted owner.
    #[must_use]
    pub fn from_bytes(bytes: Vec<u8>) -> Self {
        Self {
            bytes: Arc::new(bytes),
            origin: None,
        }
    }

    /// Reads a file, decompressing it if it is compressed.
    ///
    /// The compressed stream is consumed incrementally, so the complete
    /// compressed and decompressed buffers are never retained simultaneously.
    /// Time is `O(c + d)` for `c` compressed and `d` decompressed bytes, with
    /// `O(d)` retained memory.
    ///
    /// # Errors
    ///
    /// Returns a finding when the file cannot be read, or when decompressing it
    /// would exceed `limits`.
    pub fn open(path: impl AsRef<Path>, limits: Limits) -> Result<Self, Diagnostic> {
        let path = path.as_ref();
        let file = File::open(path)
            .map_err(|error| read_failure("input could not be opened", Some(path), error))?;

        let compressed_size = file.metadata().map(|metadata| metadata.len());
        let (compression, reader) = sniff_reader(file, path)?;

        let compressed = match compression {
            Compression::None => 0,
            Compression::Gzip | Compression::Zstd => compressed_file_size(compressed_size, path)?,
        };

        let bytes = match compression {
            Compression::None => collect_checked(
                reader,
                OutputLimit::uncompressed(limits),
                "input could not be read",
                Some(path),
            )?,
            Compression::Gzip => expand_gzip_reader(reader, compressed, limits, Some(path))?,
            Compression::Zstd => expand_zstd_reader(reader, compressed, limits, Some(path))?,
        };

        Ok(Self {
            bytes: Arc::new(bytes),
            origin: Some(path_origin(path)),
        })
    }

    /// The bytes.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        self.bytes.as_slice()
    }

    /// The number of bytes.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns true when there is nothing to read.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Where the bytes came from, for a diagnostic to name.
    ///
    /// Runs in `O(1)` time and returns a borrowed string without allocation.
    #[must_use]
    #[inline]
    pub fn origin(&self) -> Option<&str> {
        self.origin.as_deref()
    }

    /// Names where the bytes came from.
    ///
    /// Copies `origin` once into reference-counted storage in
    /// `O(origin.len())` time.
    #[must_use]
    pub fn with_origin(mut self, origin: &str) -> Self {
        self.origin = Some(Arc::from(origin));
        self
    }
}

/// A reader that replays bytes consumed while detecting compression.
///
/// The fixed-size prefix remains on the stack; no intermediate allocation is
/// required.
struct ReplayReader<R> {
    prefix: [u8; COMPRESSION_PREFIX_BYTES],
    prefix_len: usize,
    position: usize,
    inner: R,
}

impl<R> ReplayReader<R> {
    /// Creates a reader that yields `prefix[..prefix_len]` before `inner`.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[inline]
    const fn new(prefix: [u8; COMPRESSION_PREFIX_BYTES], prefix_len: usize, inner: R) -> Self {
        Self {
            prefix,
            prefix_len,
            position: 0,
            inner,
        }
    }
}

impl<R: Read> Read for ReplayReader<R> {
    /// Reads from the saved prefix first and then from the underlying stream.
    ///
    /// Each call performs `O(bytes returned)` work and allocates no memory.
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.position < self.prefix_len {
            let remaining = self.prefix_len - self.position;
            let count = remaining.min(buffer.len());

            if count == 0 {
                return Ok(0);
            }

            let end = self.position + count;
            buffer[..count].copy_from_slice(&self.prefix[self.position..end]);
            self.position = end;
            return Ok(count);
        }

        self.inner.read(buffer)
    }
}

/// The output-size policy applied while collecting bytes.
#[derive(Clone, Copy)]
enum OutputLimit {
    Uncompressed {
        maximum: u64,
    },
    Expanded {
        maximum: u64,
        compressed: u64,
        limits: Limits,
    },
}

impl OutputLimit {
    /// Creates a policy for an uncompressed stream.
    ///
    /// Only the decompressed-byte ceiling applies. Runs in `O(1)` time.
    #[inline]
    const fn uncompressed(limits: Limits) -> Self {
        Self::Uncompressed {
            maximum: limits.decompressed_bytes,
        }
    }

    /// Creates a policy for decoder output.
    ///
    /// The effective ceiling is the tighter of the absolute byte limit and the
    /// exact ratio-derived byte limit. Runs in `O(1)` time.
    #[inline]
    fn expanded(compressed: u64, limits: Limits) -> Self {
        Self::Expanded {
            maximum: maximum_expanded_bytes(compressed, limits),
            compressed,
            limits,
        }
    }

    /// Returns the largest number of output bytes this policy accepts.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[inline]
    const fn maximum(self) -> u64 {
        match self {
            Self::Uncompressed { maximum } | Self::Expanded { maximum, .. } => maximum,
        }
    }

    /// Validates a prospective output length against this policy.
    ///
    /// Returns a diagnostic before violating bytes are appended. Runs in
    /// `O(1)` time and allocates only on failure.
    #[inline]
    fn validate(self, expanded: u64) -> Result<(), Diagnostic> {
        match self {
            Self::Uncompressed { maximum } => {
                if expanded > maximum {
                    Err(Limits::exceeded("decompressed bytes", expanded))
                } else {
                    Ok(())
                }
            }
            Self::Expanded {
                compressed, limits, ..
            } => check_expansion(expanded, compressed, limits),
        }
    }
}

/// Reads the compression prefix without losing it from the stream.
///
/// Short reads and interrupted reads are handled explicitly. The returned
/// [`ReplayReader`] supplies the consumed bytes to the selected decoder before
/// continuing with the underlying reader. Runs in `O(1)` time and space.
fn sniff_reader<R: Read>(
    mut reader: R,
    path: &Path,
) -> Result<(Compression, ReplayReader<R>), Diagnostic> {
    let mut prefix = [0_u8; COMPRESSION_PREFIX_BYTES];
    let mut prefix_len = 0;

    while prefix_len < prefix.len() {
        match reader.read(&mut prefix[prefix_len..]) {
            Ok(0) => break,
            Ok(count) => {
                prefix_len += count;
            }
            Err(error) if error.kind() == ErrorKind::Interrupted => {}
            Err(error) => {
                return Err(read_failure("input could not be read", Some(path), error));
            }
        }
    }

    let compression = Compression::sniff(&prefix[..prefix_len]);
    Ok((compression, ReplayReader::new(prefix, prefix_len, reader)))
}

/// Returns the declared byte length of a compressed file.
///
/// A non-zero length is required so the compression-ratio guard cannot be
/// silently disabled by unavailable or unsuitable metadata. Runs in `O(1)`
/// time and allocates only on failure.
fn compressed_file_size(size: Result<u64, std::io::Error>, path: &Path) -> Result<u64, Diagnostic> {
    let size =
        size.map_err(|error| read_failure("input size could not be read", Some(path), error))?;

    if size == 0 {
        return Err(Diagnostic::new(Code::E1901)
            .with_message("compressed input size could not be determined")
            .with_context("path", path.display().to_string()));
    }

    Ok(size)
}

/// Expands a compressed buffer, refusing one that expands beyond the limits.
///
/// Uncompressed input is returned without copying. Compressed input is decoded
/// in `O(c + d)` time using `O(d)` output storage.
#[cfg_attr(not(test), allow(dead_code))]
fn decompress(raw: Vec<u8>, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    let compressed = match u64::try_from(raw.len()) {
        Ok(value) => value,
        Err(_) => {
            return Err(Limits::exceeded("decompressed bytes", raw.len()));
        }
    };

    match Compression::sniff(&raw) {
        Compression::None => {
            OutputLimit::uncompressed(limits).validate(compressed)?;
            Ok(raw)
        }
        Compression::Gzip => expand_gzip(&raw, compressed, limits),
        Compression::Zstd => expand_zstd(&raw, compressed, limits),
    }
}

/// Decodes a gzip slice while enforcing absolute and ratio limits.
///
/// Runs in `O(c + d)` time and uses `O(d)` output storage.
#[cfg(feature = "gzip")]
fn expand_gzip(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    expand_gzip_reader(raw, compressed, limits, None)
}

/// Reports that gzip support is not present in this build.
///
/// Runs in `O(1)` time and allocates only diagnostic context.
#[cfg(not(feature = "gzip"))]
fn expand_gzip(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("gzip"))
}

/// Decodes gzip from any reader under the configured expansion limits.
///
/// The decoder is streamed directly into the bounded collector. Runs in
/// `O(c + d)` time and retains only the `O(d)` output buffer.
#[cfg(feature = "gzip")]
fn expand_gzip_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    let decoder = flate2::read::GzDecoder::new(reader);

    collect_checked(
        decoder,
        OutputLimit::expanded(compressed, limits),
        "gzip stream could not be decoded",
        path,
    )
}

/// Reports that gzip support is not present in this build.
///
/// Runs in `O(1)` time and allocates only diagnostic context.
#[cfg(not(feature = "gzip"))]
fn expand_gzip_reader<R: Read>(
    _: R,
    _: u64,
    _: Limits,
    _: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("gzip"))
}

/// Decodes a Zstandard slice while enforcing absolute and ratio limits.
///
/// Runs in `O(c + d)` time and uses `O(d)` output storage.
#[cfg(feature = "zstd")]
fn expand_zstd(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    expand_zstd_reader(raw, compressed, limits, None)
}

/// Reports that Zstandard support is not present in this build.
///
/// Runs in `O(1)` time and allocates only diagnostic context.
#[cfg(not(feature = "zstd"))]
fn expand_zstd(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("zstd"))
}

/// Decodes Zstandard from any reader under the configured expansion limits.
///
/// The decoder is streamed directly into the bounded collector. Runs in
/// `O(c + d)` time and retains only the `O(d)` output buffer.
#[cfg(feature = "zstd")]
fn expand_zstd_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    let decoder = zstd::stream::read::Decoder::new(reader)
        .map_err(|error| read_failure("zstd stream could not be decoded", path, error))?;

    collect_checked(
        decoder,
        OutputLimit::expanded(compressed, limits),
        "zstd stream could not be decoded",
        path,
    )
}

/// Reports that Zstandard support is not present in this build.
///
/// Runs in `O(1)` time and allocates only diagnostic context.
#[cfg(not(feature = "zstd"))]
fn expand_zstd_reader<R: Read>(
    _: R,
    _: u64,
    _: Limits,
    _: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("zstd"))
}

/// Collects bytes while validating every prospective output length.
///
/// Validation occurs before appending a chunk, so the retained length never
/// exceeds `limit`. Reads are restricted to the remaining allowance plus one
/// byte, preventing a decoder from performing an entire large output read after
/// the limit has already been reached.
///
/// Capacity grows geometrically and is capped by the effective output ceiling.
/// Collection takes amortized `O(n)` time and `O(n)` output storage.
fn collect_checked<R: Read>(
    mut reader: R,
    limit: OutputLimit,
    failure_message: &'static str,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    let mut output = Vec::new();
    let mut buffer = [0_u8; READ_BUFFER_BYTES];

    loop {
        let current = match u64::try_from(output.len()) {
            Ok(value) => value,
            Err(_) => {
                return Err(Limits::exceeded("decompressed bytes", output.len()));
            }
        };

        let remaining = limit.maximum().saturating_sub(current);
        let request = match usize::try_from(remaining) {
            Ok(value) if value < buffer.len() => value + 1,
            Ok(_) | Err(_) => buffer.len(),
        };

        let count = match reader.read(&mut buffer[..request]) {
            Ok(0) => break,
            Ok(count) => count,
            Err(error) if error.kind() == ErrorKind::Interrupted => {
                continue;
            }
            Err(error) => {
                return Err(read_failure(failure_message, path, error));
            }
        };

        let additional = match u64::try_from(count) {
            Ok(value) => value,
            Err(_) => {
                return Err(Limits::exceeded("decompressed bytes", count));
            }
        };

        let expanded = match current.checked_add(additional) {
            Some(value) => value,
            None => {
                return Err(Limits::exceeded("decompressed bytes", "more than u64::MAX"));
            }
        };

        limit.validate(expanded)?;
        reserve_output(&mut output, count, limit.maximum(), path)?;
        output.extend_from_slice(&buffer[..count]);
    }

    Ok(output)
}

/// Reserves enough output capacity for one append without unbounded overgrowth.
///
/// Capacity grows geometrically for amortized linear collection, but never
/// intentionally beyond the effective output ceiling. Allocation and capacity
/// failures become diagnostics rather than indexing or capacity panics.
fn reserve_output(
    output: &mut Vec<u8>,
    additional: usize,
    maximum: u64,
    path: Option<&Path>,
) -> Result<(), Diagnostic> {
    let required = match output.len().checked_add(additional) {
        Some(value) => value,
        None => {
            return Err(Limits::exceeded(
                "decompressed bytes",
                "more than usize::MAX",
            ));
        }
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

    let grown = if output.capacity() == 0 {
        required
    } else {
        output.capacity().saturating_mul(2).max(required)
    };
    let target = grown.min(maximum_capacity);
    let extra = match target.checked_sub(output.len()) {
        Some(value) => value,
        None => {
            return Err(Limits::exceeded("decompressed bytes", required));
        }
    };

    output
        .try_reserve_exact(extra)
        .map_err(|error| allocation_failure(path, error))
}

/// Builds a diagnostic for a file or decoder read failure.
///
/// The path is included when available. Formatting takes
/// `O(path length + reason length)` time and allocates only diagnostic context.
fn read_failure(message: &'static str, path: Option<&Path>, error: std::io::Error) -> Diagnostic {
    let diagnostic = Diagnostic::new(Code::E1901).with_message(message);

    let diagnostic = match path {
        Some(path) => diagnostic.with_context("path", path.display().to_string()),
        None => diagnostic,
    };

    diagnostic.with_context("reason", error.to_string())
}

/// Builds a diagnostic for an output-buffer allocation failure.
///
/// The path is included when available. This allocates only diagnostic context.
fn allocation_failure(path: Option<&Path>, error: TryReserveError) -> Diagnostic {
    let diagnostic =
        Diagnostic::new(Code::E1901).with_message("input buffer could not be allocated");

    let diagnostic = match path {
        Some(path) => diagnostic.with_context("path", path.display().to_string()),
        None => diagnostic,
    };

    diagnostic.with_context("reason", error.to_string())
}

/// Converts a path into shared diagnostic-origin storage.
///
/// Valid UTF-8 paths are copied directly; non-UTF-8 paths use the same lossy
/// representation as `Path::display`. Runs in `O(path.len())` time.
fn path_origin(path: &Path) -> Arc<str> {
    match path.to_str() {
        Some(origin) => Arc::from(origin),
        None => Arc::from(path.to_string_lossy().as_ref()),
    }
}

/// Returns the effective maximum decoder-output length.
///
/// The result is the minimum of the absolute decompressed-byte ceiling and the
/// exact ratio-derived ceiling. Runs in `O(1)` time and allocates no memory.
fn maximum_expanded_bytes(compressed: u64, limits: Limits) -> u64 {
    limits
        .decompressed_bytes
        .min(ratio_byte_limit(compressed, limits.compression_ratio))
}

/// Returns the largest output size accepted by a compression ratio.
///
/// Multiplication overflow means the mathematical limit exceeds every
/// representable `u64` output size. A zero compressed size disables this
/// calculation defensively; file-backed compressed input rejects such metadata
/// before decoding. Runs in `O(1)` time.
fn ratio_byte_limit(compressed: u64, ratio: u64) -> u64 {
    if compressed == 0 {
        return u64::MAX;
    }

    match compressed.checked_mul(ratio) {
        Some(value) => value,
        None => u64::MAX,
    }
}

/// Reports that a compression implementation is absent from this build.
///
/// Runs in `O(1)` time and allocates only diagnostic context.
#[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
fn unsupported(container: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message("input is compressed with a container this build does not include")
        .with_context("container", container)
}

/// Checks absolute decompressed size and exact expansion ratio.
///
/// The ratio comparison uses multiplication rather than integer division, so a
/// fractional overrun cannot be rounded down and accepted. Runs in `O(1)` time
/// and allocates only when constructing a failure diagnostic.
#[cfg_attr(
    not(any(feature = "gzip", feature = "zstd")),
    expect(dead_code, reason = "only a compression feature reaches this")
)]
fn check_expansion(expanded: u64, compressed: u64, limits: Limits) -> Result<(), Diagnostic> {
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

/// Returns `ceil(numerator / denominator)`.
///
/// Runs in `O(1)` time and allocates no memory. A zero denominator produces
/// `u64::MAX` defensively, although file-backed compressed inputs reject that
/// condition before decoding.
fn ceiling_ratio(numerator: u64, denominator: u64) -> u64 {
    if denominator == 0 {
        return u64::MAX;
    }

    let quotient = numerator / denominator;

    if numerator % denominator == 0 {
        quotient
    } else {
        match quotient.checked_add(1) {
            Some(value) => value,
            None => u64::MAX,
        }
    }
}

#[cfg(test)]
#[path = "input_tests.rs"]
mod tests;
