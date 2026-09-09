//! Bounded, checksummed temporary records owned by one execution graph.

use super::TempStoragePolicy;
use std::fmt;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

const FILE_MAGIC: [u8; 8] = *b"PDBXSP01";
pub(crate) const FILE_HEADER_BYTES: u64 = FILE_MAGIC.len() as u64;
const RECORD_HEADER_BYTES: usize = 16;
pub(crate) const RECORD_HEADER_DISK_BYTES: u64 = 16;
const RECORD_HEADER_REWIND: i64 = -16;
const HASH_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const HASH_PRIME: u64 = 0x0000_0100_0000_01b3;

#[derive(Debug)]
pub(crate) struct DiskAccount {
    root: PathBuf,
    maximum: u64,
    used: AtomicU64,
    spilled: AtomicU64,
}

impl DiskAccount {
    pub(crate) fn from_policy(policy: &TempStoragePolicy) -> Option<Arc<Self>> {
        match policy {
            TempStoragePolicy::Disabled => None,
            TempStoragePolicy::Directory { root, max_bytes } => Some(Arc::new(Self {
                root: root.clone(),
                maximum: *max_bytes,
                used: AtomicU64::new(0),
                spilled: AtomicU64::new(0),
            })),
        }
    }

    fn try_charge(&self, bytes: u64) -> Result<(), SpillError> {
        let mut current = self.used.load(Ordering::Acquire);
        loop {
            let Some(next) = current.checked_add(bytes) else {
                return Err(self.exhausted(bytes, current));
            };
            if next > self.maximum {
                return Err(self.exhausted(bytes, current));
            }
            match self.used.compare_exchange_weak(
                current,
                next,
                Ordering::AcqRel,
                Ordering::Acquire,
            ) {
                Ok(_) => return Ok(()),
                Err(observed) => current = observed,
            }
        }
    }

    fn exhausted(&self, requested: u64, used: u64) -> SpillError {
        SpillError::BudgetExceeded {
            requested,
            available: self.maximum.saturating_sub(used),
        }
    }

    fn release(&self, bytes: u64) {
        self.used.fetch_sub(bytes, Ordering::AcqRel);
    }

    pub(crate) fn used_bytes(&self) -> u64 {
        self.used.load(Ordering::Acquire)
    }

    pub(crate) fn spilled_bytes(&self) -> u64 {
        self.spilled.load(Ordering::Acquire)
    }
}

/// Failure to create, write, or validate bounded temporary storage.
#[derive(Debug)]
pub enum SpillError {
    /// The execution context does not permit temporary disk storage.
    Disabled,
    /// A key must be one portable file-name component.
    InvalidKey,
    /// The configured disk ceiling cannot retain the requested bytes.
    BudgetExceeded {
        /// Additional bytes requested.
        requested: u64,
        /// Unreserved bytes remaining.
        available: u64,
    },
    /// A record length cannot be represented on this platform.
    RecordTooLarge(u64),
    /// The next record exceeds the caller's explicit memory ceiling.
    RecordExceedsLimit {
        /// Bytes declared by the record.
        record: u64,
        /// Bytes the caller is prepared to retain.
        limit: usize,
    },
    /// The file is not a supported pdbiox spill artifact.
    CorruptHeader,
    /// A record was truncated or failed its checksum.
    CorruptRecord,
    /// The operating system rejected an I/O operation.
    Io(io::Error),
}

impl fmt::Display for SpillError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Disabled => formatter.write_str("temporary disk storage is disabled"),
            Self::InvalidKey => formatter.write_str("spill key must be one portable file name"),
            Self::BudgetExceeded {
                requested,
                available,
            } => write!(
                formatter,
                "spill requested {requested} bytes with {available} bytes available"
            ),
            Self::RecordTooLarge(bytes) => {
                write!(
                    formatter,
                    "spill record of {bytes} bytes is not addressable"
                )
            }
            Self::RecordExceedsLimit { record, limit } => {
                write!(
                    formatter,
                    "spill record has {record} bytes but limit is {limit}"
                )
            }
            Self::CorruptHeader => formatter.write_str("spill file header is invalid"),
            Self::CorruptRecord => formatter.write_str("spill record is truncated or corrupt"),
            Self::Io(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for SpillError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for SpillError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// Append-only temporary file charged to one execution context.
#[derive(Debug)]
pub struct SpillFile {
    file: Option<File>,
    path: Option<PathBuf>,
    account: Arc<DiskAccount>,
    charged: u64,
    records: u64,
}

impl SpillFile {
    pub(crate) fn create(account: Arc<DiskAccount>, key: &str) -> Result<Self, SpillError> {
        if !portable_key(key) {
            return Err(SpillError::InvalidKey);
        }
        fs::create_dir_all(&account.root)?;
        account.try_charge(FILE_HEADER_BYTES)?;
        let path = account.root.join(format!("{key}.pdbiox-spill"));
        let opened = OpenOptions::new().write(true).create_new(true).open(&path);
        let mut file = match opened {
            Ok(file) => file,
            Err(error) => {
                account.release(FILE_HEADER_BYTES);
                return Err(SpillError::Io(error));
            }
        };
        if let Err(error) = file.write_all(&FILE_MAGIC) {
            account.release(FILE_HEADER_BYTES);
            let _ignored = fs::remove_file(&path);
            return Err(SpillError::Io(error));
        }
        account
            .spilled
            .fetch_add(FILE_HEADER_BYTES, Ordering::AcqRel);
        Ok(Self {
            file: Some(file),
            path: Some(path),
            account,
            charged: FILE_HEADER_BYTES,
            records: 0,
        })
    }

    /// Appends one checksummed record without allocating.
    ///
    /// # Errors
    ///
    /// Returns an error without retaining a partial record or disk charge.
    pub fn append(&mut self, payload: &[u8]) -> Result<(), SpillError> {
        let payload_bytes =
            u64::try_from(payload.len()).map_err(|_error| SpillError::RecordTooLarge(u64::MAX))?;
        let Some(frame_bytes) = payload_bytes.checked_add(RECORD_HEADER_DISK_BYTES) else {
            return Err(SpillError::RecordTooLarge(payload_bytes));
        };
        self.account.try_charge(frame_bytes)?;
        let Some(file) = &mut self.file else {
            self.account.release(frame_bytes);
            return Err(SpillError::CorruptHeader);
        };
        let original_len = self.charged;
        let mut header = [0_u8; RECORD_HEADER_BYTES];
        header[..8].copy_from_slice(&payload_bytes.to_le_bytes());
        header[8..].copy_from_slice(&checksum(payload).to_le_bytes());
        if let Err(error) = file
            .write_all(&header)
            .and_then(|()| file.write_all(payload))
        {
            let _ignored = file.set_len(original_len);
            let _ignored = file.seek(SeekFrom::End(0));
            self.account.release(frame_bytes);
            return Err(SpillError::Io(error));
        }
        self.charged += frame_bytes;
        self.records += 1;
        self.account
            .spilled
            .fetch_add(frame_bytes, Ordering::AcqRel);
        Ok(())
    }

    /// Flushes the writer and returns an immutable temporary artifact.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if buffered bytes cannot reach the filesystem.
    pub fn finish(mut self) -> Result<SpillArtifact, SpillError> {
        let Some(mut file) = self.file.take() else {
            return Err(SpillError::CorruptHeader);
        };
        file.flush()?;
        file.sync_data()?;
        drop(file);
        let Some(path) = self.path.take() else {
            return Err(SpillError::CorruptHeader);
        };
        let artifact = SpillArtifact {
            path,
            account: Arc::clone(&self.account),
            charged: self.charged,
            records: self.records,
        };
        self.charged = 0;
        Ok(artifact)
    }
}

impl Drop for SpillFile {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            let _ignored = fs::remove_file(path);
        }
        self.account.release(self.charged);
    }
}

/// Immutable temporary records deleted when the artifact is dropped.
#[derive(Debug)]
pub struct SpillArtifact {
    path: PathBuf,
    account: Arc<DiskAccount>,
    charged: u64,
    records: u64,
}

impl SpillArtifact {
    /// Opens an independent sequential reader and validates the file header.
    ///
    /// # Errors
    ///
    /// Returns an error if the artifact is absent or has an invalid header.
    pub fn reader(&self) -> Result<SpillReader, SpillError> {
        SpillReader::open(&self.path)
    }

    /// Returns the physical bytes charged to the execution context.
    #[must_use]
    pub const fn bytes(&self) -> u64 {
        self.charged
    }

    /// Returns the number of complete records written.
    #[must_use]
    pub const fn records(&self) -> u64 {
        self.records
    }

    /// Returns the artifact path for diagnostics and external disk monitoring.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for SpillArtifact {
    fn drop(&mut self) {
        let _ignored = fs::remove_file(&self.path);
        self.account.release(self.charged);
    }
}

/// Sequential reader that reuses caller-owned payload capacity.
#[derive(Debug)]
pub struct SpillReader {
    file: File,
}

impl SpillReader {
    fn open(path: &Path) -> Result<Self, SpillError> {
        let mut file = File::open(path)?;
        let mut magic = [0_u8; FILE_MAGIC.len()];
        file.read_exact(&mut magic)
            .map_err(|_error| SpillError::CorruptHeader)?;
        if magic != FILE_MAGIC {
            return Err(SpillError::CorruptHeader);
        }
        Ok(Self { file })
    }

    /// Reads one record into a reusable buffer. `None` marks clean end-of-file.
    ///
    /// # Errors
    ///
    /// Returns [`SpillError::CorruptRecord`] for truncation or checksum drift.
    pub fn read_next(
        &mut self,
        payload: &mut Vec<u8>,
        max_bytes: usize,
    ) -> Result<Option<usize>, SpillError> {
        let mut header = [0_u8; RECORD_HEADER_BYTES];
        let mut first = [0_u8; 1];
        if self.file.read(&mut first)? == 0 {
            return Ok(None);
        }
        header[0] = first[0];
        self.file
            .read_exact(&mut header[1..])
            .map_err(|_error| SpillError::CorruptRecord)?;
        let length_u64 = u64::from_le_bytes(copy_u64(&header[..8]));
        let expected = u64::from_le_bytes(copy_u64(&header[8..]));
        let payload_offset = self.file.stream_position()?;
        let file_bytes = self.file.metadata()?.len();
        if length_u64 > file_bytes.saturating_sub(payload_offset) {
            return Err(SpillError::CorruptRecord);
        }
        let limit_u64 = match u64::try_from(max_bytes) {
            Ok(limit) => limit,
            Err(_error) => u64::MAX,
        };
        if length_u64 > limit_u64 {
            self.file.seek(SeekFrom::Current(RECORD_HEADER_REWIND))?;
            return Err(SpillError::RecordExceedsLimit {
                record: length_u64,
                limit: max_bytes,
            });
        }
        let length =
            usize::try_from(length_u64).map_err(|_error| SpillError::RecordTooLarge(length_u64))?;
        payload.resize(length, 0);
        self.file
            .read_exact(payload)
            .map_err(|_error| SpillError::CorruptRecord)?;
        if checksum(payload) != expected {
            return Err(SpillError::CorruptRecord);
        }
        Ok(Some(length))
    }
}

fn portable_key(key: &str) -> bool {
    !key.is_empty()
        && key.len() <= 96
        && key
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_'))
}

pub(crate) fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(HASH_OFFSET, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(HASH_PRIME)
    })
}

fn copy_u64(bytes: &[u8]) -> [u8; 8] {
    let mut value = [0_u8; 8];
    value.copy_from_slice(bytes);
    value
}

#[cfg(test)]
#[path = "spill_tests.rs"]
mod tests;
