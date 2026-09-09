//! Owned read-only snapshots and explicitly unsafe file-backed mappings.

use memmap2::{Advice, Mmap, MmapMut, MmapOptions};
use std::fs::File;
use std::io::{self, Read, Seek, Write};
use std::ops::Deref;

#[derive(Debug)]
enum Mapping {
    Empty,
    Snapshot(Mmap),
    FileBacked { mapping: Mmap, _file: File },
}

/// Read-only bytes backed by an owned snapshot or an unchecked file mapping.
///
/// [`Self::new`] and [`Self::snapshot`] are safe because they copy the current
/// contents into private anonymous memory. Zero-copy access to an external file
/// is available only through [`Self::map_file_unchecked`], whose caller must
/// uphold the backing-file stability contract required by `memmap2`.
#[derive(Debug)]
pub struct MappedFile {
    mapping: Mapping,
}

impl MappedFile {
    /// Streams bytes into a private temporary file and maps the finished snapshot.
    ///
    /// The temporary object has no caller-visible path and its only retained
    /// descriptor is owned by the returned mapping. This keeps resident memory
    /// proportional to the copy buffer while preserving immutable byte-slice
    /// semantics for parsers.
    ///
    /// # Errors
    ///
    /// Returns an I/O error when the temporary file cannot be created, filled,
    /// rewound, or mapped.
    pub fn private_snapshot_from_reader(mut reader: impl Read) -> io::Result<Self> {
        let mut file = tempfile::tempfile()?;
        let mut buffer = vec![0_u8; 64 * 1024].into_boxed_slice();
        loop {
            let count = match reader.read(&mut buffer) {
                Ok(0) => break,
                Ok(count) => count,
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            };
            file.write_all(&buffer[..count])?;
        }
        file.flush()?;
        file.rewind()?;
        // SAFETY: `tempfile()` creates a private object without a caller-visible
        // path. No descriptor escaped before this mapping, and the mapping
        // retains its sole live descriptor for its complete lifetime.
        unsafe { Self::map_file_unchecked(&file) }
    }

    /// Creates a private, read-only snapshot of the complete current file.
    ///
    /// This compatibility constructor is equivalent to [`Self::snapshot`]. Use
    /// the named constructor when the distinction from an unchecked zero-copy
    /// file mapping should be explicit at the call site.
    ///
    /// # Errors
    ///
    /// Returns the same operating-system errors as [`Self::snapshot`].
    pub fn new(file: &File) -> io::Result<Self> {
        Self::snapshot(file)
    }

    /// Copies the complete current file into private, read-only anonymous pages.
    ///
    /// A concurrent writer can make the read fail or produce bytes assembled
    /// from different writes, but it cannot mutate or invalidate the returned
    /// snapshot.
    ///
    /// # Errors
    ///
    /// Returns an operating-system error when the length cannot be represented,
    /// the complete snapshot cannot be read, or anonymous mapping fails.
    pub fn snapshot(file: &File) -> io::Result<Self> {
        let length = file_length(file)?;
        if length == 0 {
            return Ok(Self {
                mapping: Mapping::Empty,
            });
        }

        let mut snapshot = MmapMut::map_anon(length)?;
        read_at(file, &mut snapshot)?;
        let mapping = snapshot.make_read_only()?;
        Ok(Self {
            mapping: Mapping::Snapshot(mapping),
        })
    }

    /// Maps a file directly without copying it.
    ///
    /// The returned value retains a cloned descriptor, so its lifetime does not
    /// borrow `file`. Descriptor ownership alone does not make the mapping safe:
    /// other handles and processes can still mutate the backing object.
    ///
    /// # Safety
    ///
    /// The caller must guarantee that the mapped byte range is not modified,
    /// truncated, or otherwise invalidated for the complete lifetime of the
    /// returned value. The guarantee must cover every handle and process able to
    /// mutate the backing object. Read-only permissions on this descriptor and
    /// retaining it inside the returned value do not establish that invariant.
    ///
    /// # Errors
    ///
    /// Returns an operating-system error when the descriptor cannot be cloned,
    /// the length cannot be represented, or the mapping cannot be created.
    pub unsafe fn map_file_unchecked(file: &File) -> io::Result<Self> {
        let length = file_length(file)?;
        if length == 0 {
            return Ok(Self {
                mapping: Mapping::Empty,
            });
        }

        let retained = file.try_clone()?;
        // SAFETY: the caller accepts the complete file-stability contract
        // documented on this constructor for the returned mapping's lifetime.
        let mapping = unsafe { MmapOptions::new().len(length).map(&retained)? };
        Ok(Self {
            mapping: Mapping::FileBacked {
                mapping,
                _file: retained,
            },
        })
    }

    /// Tells the kernel this mapping will be read front to back, once.
    ///
    /// Every parser above this crate walks its input strictly forwards and
    /// never returns to it. Saying so lets the kernel read ahead further and
    /// stop retaining pages that have been passed, which is the difference
    /// between a hundred-gigabyte scan that fits in page cache and one that
    /// evicts everything else on the machine.
    ///
    /// Advice is a hint, and a platform that declines it has not failed to
    /// read the file. Callers are expected to ignore a refusal.
    ///
    /// # Errors
    ///
    /// Returns the operating system's error when it refuses the advice.
    pub fn advise_sequential(&self) -> io::Result<()> {
        match &self.mapping {
            Mapping::Empty => Ok(()),
            Mapping::Snapshot(mapping) | Mapping::FileBacked { mapping, .. } => {
                mapping.advise(Advice::Sequential)
            }
        }
    }

    /// The mapped bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        match &self.mapping {
            Mapping::Empty => &[],
            Mapping::Snapshot(mapping) | Mapping::FileBacked { mapping, .. } => mapping,
        }
    }
}

impl AsRef<[u8]> for MappedFile {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Deref for MappedFile {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.as_bytes()
    }
}

fn file_length(file: &File) -> io::Result<usize> {
    usize::try_from(file.metadata()?.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "file length exceeds the platform address space",
        )
    })
}

#[cfg(unix)]
fn read_at(file: &File, destination: &mut [u8]) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    read_chunks(destination, |buffer, offset| file.read_at(buffer, offset))
}

#[cfg(windows)]
fn read_at(file: &File, destination: &mut [u8]) -> io::Result<()> {
    use std::os::windows::fs::FileExt;

    read_chunks(destination, |buffer, offset| file.seek_read(buffer, offset))
}

#[cfg(not(any(unix, windows)))]
fn read_at(file: &File, destination: &mut [u8]) -> io::Result<()> {
    use std::io::{Read, Seek};

    let mut retained = file.try_clone()?;
    retained.rewind()?;
    retained.read_exact(destination)
}

#[cfg(any(unix, windows))]
fn read_chunks(
    mut destination: &mut [u8],
    mut read: impl FnMut(&mut [u8], u64) -> io::Result<usize>,
) -> io::Result<()> {
    let mut offset = 0_u64;
    while !destination.is_empty() {
        let count = read(destination, offset)?;
        if count == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "file changed while its snapshot was being read",
            ));
        }
        offset = offset
            .checked_add(
                u64::try_from(count).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "read size exceeds u64")
                })?,
            )
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "read offset overflow"))?;
        destination = &mut destination[count..];
    }
    Ok(())
}

#[cfg(test)]
#[path = "mapped_tests.rs"]
mod tests;
