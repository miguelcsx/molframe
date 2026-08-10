//! Owned read-only snapshots and explicitly unsafe file-backed mappings.

use memmap2::{Mmap, MmapMut, MmapOptions};
use std::fs::File;
use std::io;
use std::ops::Deref;

#[derive(Debug)]
enum Mapping {
    Empty,
    Snapshot(Mmap),
    FileBacked { mapping: Mmap, _file: File },
}

/// Read-only bytes backed by an owned snapshot or an explicitly unsafe file map.
///
/// [`Self::new`] is the safe default: it copies the current file contents into
/// anonymous memory before publishing a read-only mapping. Later changes to the
/// source file cannot invalidate or alter that snapshot.
#[derive(Debug)]
pub struct MappedFile {
    mapping: Mapping,
}

impl MappedFile {
    /// Copies the complete current file into an owned, read-only anonymous map.
    ///
    /// This is a snapshot operation, not a zero-copy file mapping. A concurrent
    /// writer can make the read fail or produce bytes from different writes, but
    /// cannot cause undefined behaviour or mutate the returned snapshot.
    ///
    /// # Errors
    ///
    /// Returns an operating-system error when the length cannot be represented,
    /// the complete snapshot cannot be read, or anonymous mapping fails.
    pub fn new(file: &File) -> io::Result<Self> {
        let length = usize::try_from(file.metadata()?.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "file length exceeds the platform address space",
            )
        })?;
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
    /// This operation is separated from [`Self::new`] because Rust cannot stop
    /// another handle, process, or filesystem actor from mutating a file-backed
    /// mapping.
    ///
    /// # Safety
    ///
    /// The caller must guarantee that the mapped byte range is not modified,
    /// truncated, replaced through the same inode, or otherwise invalidated for
    /// the complete lifetime of the returned value. That guarantee must cover
    /// every handle and process able to mutate the backing object. Read-only
    /// mapping permissions and retaining this descriptor do not establish that
    /// invariant by themselves.
    ///
    /// # Errors
    ///
    /// Returns the operating-system error when the descriptor cannot be cloned
    /// or the mapping cannot be created.
    pub unsafe fn map_file_unchecked(file: &File) -> io::Result<Self> {
        let retained = file.try_clone()?;
        // SAFETY: the caller accepts the complete file-stability contract
        // documented on this unsafe constructor for the returned lifetime.
        let mapping = unsafe { MmapOptions::new().map(file)? };
        Ok(Self {
            mapping: Mapping::FileBacked {
                mapping,
                _file: retained,
            },
        })
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
