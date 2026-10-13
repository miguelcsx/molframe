//! Frame-offset indexes for formats whose frames vary in length.
//!
//! Compressed XTC frames and TRR frames with optional velocity or force blocks
//! cannot be located by arithmetic, so random access needs the byte offset of
//! every frame. One sequential header scan produces it. The index is cheap to
//! share between readers (an `Arc` of offsets) and serialises to a small
//! checksummed blob, so a second process, or a second run, can skip the scan.
//! Deserialising validates every declared size against the bytes actually
//! present before allocating, and checks the stamp of the source file the
//! index was built from, so a stale or damaged cache is refused, never trusted.

use std::io::ErrorKind;
use std::path::Path;
use std::sync::Arc;
use std::time::UNIX_EPOCH;

use crate::TrajectoryError;

const MAGIC: [u8; 4] = *b"MFIX";
const VERSION: u32 = 1;
const HEADER_BYTES: usize = 48;
const CHECKSUM_BYTES: usize = 8;
const OFFSET_BYTES: usize = size_of::<u64>();

/// Identity of the file an index describes: its length and modification time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FileStamp {
    /// Length in bytes.
    pub len: u64,
    /// Modification time in nanoseconds since the Unix epoch; zero when unknown.
    pub modified_ns: u64,
}

impl FileStamp {
    /// Reads the stamp of `path`.
    ///
    /// # Errors
    ///
    /// Returns [`FrameIndexError::Io`] when the file cannot be inspected.
    pub fn of(path: &Path) -> Result<Self, FrameIndexError> {
        let metadata =
            std::fs::metadata(path).map_err(|error| FrameIndexError::Io(error.kind()))?;
        Ok(Self::from_metadata(&metadata))
    }

    pub(crate) fn from_metadata(metadata: &std::fs::Metadata) -> Self {
        let modified = metadata
            .modified()
            .ok()
            .and_then(|time| time.duration_since(UNIX_EPOCH).ok())
            .and_then(|elapsed| u64::try_from(elapsed.as_nanos()).ok());
        let Some(modified_ns) = modified else {
            return Self {
                len: metadata.len(),
                modified_ns: 0,
            };
        };
        Self {
            len: metadata.len(),
            modified_ns,
        }
    }
}

/// An invalid, stale or unreadable frame index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum FrameIndexError {
    /// The index file could not be read or written.
    #[error("frame index I/O failed: {0:?}")]
    Io(ErrorKind),
    /// The bytes are not a frame index, are truncated, or fail validation.
    #[error("frame index is malformed")]
    Malformed,
    /// The index describes a different file, format or atom count.
    #[error("frame index does not describe this trajectory")]
    Mismatch,
}

impl From<FrameIndexError> for TrajectoryError {
    fn from(error: FrameIndexError) -> Self {
        match error {
            FrameIndexError::Io(kind) => Self::SourceIo {
                format: "frame index",
                kind,
            },
            FrameIndexError::Malformed | FrameIndexError::Mismatch => Self::InvalidSource {
                format: "frame index",
            },
        }
    }
}

/// Byte offset of every frame in one trajectory file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameIndex {
    format: [u8; 4],
    atoms: u64,
    stamp: FileStamp,
    offsets: Arc<[u64]>,
}

impl FrameIndex {
    /// Format tag of XTC indexes.
    pub const XTC: [u8; 4] = *b"xtc ";
    /// Format tag of TRR indexes.
    pub const TRR: [u8; 4] = *b"trr ";

    /// Builds an index from scanned offsets.
    ///
    /// # Errors
    ///
    /// Returns [`FrameIndexError::Malformed`] unless the offsets start at zero,
    /// strictly increase and stay inside the file.
    pub fn new(
        format: [u8; 4],
        atoms: usize,
        stamp: FileStamp,
        offsets: Vec<u64>,
    ) -> Result<Self, FrameIndexError> {
        let ordered = offsets.first().is_none_or(|first| *first == 0)
            && offsets.windows(2).all(|pair| pair[0] < pair[1])
            && offsets.last().is_none_or(|last| *last < stamp.len);
        if !ordered {
            return Err(FrameIndexError::Malformed);
        }
        Ok(Self {
            format,
            atoms: u64::try_from(atoms).map_err(|_| FrameIndexError::Malformed)?,
            stamp,
            offsets: offsets.into(),
        })
    }

    /// Frames indexed.
    #[must_use]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Whether no frame is indexed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Byte offset of each frame, in file order.
    #[must_use]
    pub fn offsets(&self) -> &[u64] {
        &self.offsets
    }

    pub(crate) fn shared(&self) -> Arc<[u64]> {
        Arc::clone(&self.offsets)
    }

    /// Stamp of the file the index was built from.
    #[must_use]
    pub const fn stamp(&self) -> FileStamp {
        self.stamp
    }

    /// Whether the index was built for this format, atom count and file.
    #[must_use]
    pub fn describes(&self, format: [u8; 4], atoms: usize, stamp: FileStamp) -> bool {
        self.format == format
            && u64::try_from(atoms).is_ok_and(|atoms| atoms == self.atoms)
            && self.stamp == stamp
    }

    /// Serialises the index with a trailing checksum.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes =
            Vec::with_capacity(HEADER_BYTES + self.offsets.len() * OFFSET_BYTES + CHECKSUM_BYTES);
        bytes.extend_from_slice(&MAGIC);
        bytes.extend_from_slice(&VERSION.to_le_bytes());
        bytes.extend_from_slice(&self.format);
        bytes.extend_from_slice(&[0; 4]);
        bytes.extend_from_slice(&self.atoms.to_le_bytes());
        bytes.extend_from_slice(&self.stamp.len.to_le_bytes());
        bytes.extend_from_slice(&self.stamp.modified_ns.to_le_bytes());
        bytes.extend_from_slice(&(self.offsets.len() as u64).to_le_bytes());
        for offset in self.offsets.iter() {
            bytes.extend_from_slice(&offset.to_le_bytes());
        }
        let sum = checksum(&bytes);
        bytes.extend_from_slice(&sum.to_le_bytes());
        bytes
    }

    /// Parses and validates a serialised index.
    ///
    /// The declared frame count is checked against the bytes present before
    /// anything is allocated, so a forged count cannot reserve memory.
    ///
    /// # Errors
    ///
    /// Returns [`FrameIndexError::Malformed`] for a bad magic, version,
    /// length, checksum or offset order.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, FrameIndexError> {
        let body_len = bytes
            .len()
            .checked_sub(CHECKSUM_BYTES)
            .filter(|len| *len >= HEADER_BYTES)
            .ok_or(FrameIndexError::Malformed)?;
        let (body, tail) = bytes.split_at(body_len);
        let stored = u64::from_le_bytes(tail.try_into().map_err(|_| FrameIndexError::Malformed)?);
        if stored != checksum(body) || body[..4] != MAGIC || word32(body, 4) != VERSION {
            return Err(FrameIndexError::Malformed);
        }
        let count = usize::try_from(word64(body, 40)).map_err(|_| FrameIndexError::Malformed)?;
        if count.checked_mul(OFFSET_BYTES) != Some(body_len - HEADER_BYTES) {
            return Err(FrameIndexError::Malformed);
        }
        let offsets = body[HEADER_BYTES..]
            .as_chunks::<OFFSET_BYTES>()
            .0
            .iter()
            .map(|chunk| u64::from_le_bytes(*chunk))
            .collect();
        let format = [body[8], body[9], body[10], body[11]];
        let atoms = usize::try_from(word64(body, 16)).map_err(|_| FrameIndexError::Malformed)?;
        let stamp = FileStamp {
            len: word64(body, 24),
            modified_ns: word64(body, 32),
        };
        Self::new(format, atoms, stamp, offsets)
    }

    /// Writes the serialised index to `path`.
    ///
    /// # Errors
    ///
    /// Returns [`FrameIndexError::Io`] when the file cannot be written.
    pub fn write(&self, path: &Path) -> Result<(), FrameIndexError> {
        std::fs::write(path, self.to_bytes()).map_err(|error| FrameIndexError::Io(error.kind()))
    }

    /// Reads and validates an index file.
    ///
    /// # Errors
    ///
    /// Returns an I/O error, or the validation errors of [`Self::from_bytes`].
    pub fn read(path: &Path) -> Result<Self, FrameIndexError> {
        let bytes = std::fs::read(path).map_err(|error| FrameIndexError::Io(error.kind()))?;
        Self::from_bytes(&bytes)
    }
}

fn word32(bytes: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
}

fn word64(bytes: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(array8(&bytes[at..at + 8]))
}

fn array8(chunk: &[u8]) -> [u8; 8] {
    let mut out = [0; 8];
    out.copy_from_slice(&chunk[..8]);
    out
}

/// FNV-1a over 64 bits: detects truncation and bit rot, not tampering.
fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

#[cfg(test)]
#[path = "frame_index_tests.rs"]
mod tests;
