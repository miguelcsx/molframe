//! Ownership wrapper around a read-only file mapping.

use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::io;
use std::ops::Deref;

/// A read-only mapping that keeps its source file open for its full lifetime.
///
/// No mutable view or source handle escapes this type. Keeping the descriptor
/// alive also avoids mappings whose backing object is closed by this process
/// while a parser still owns the bytes.
#[derive(Debug)]
pub struct MappedFile {
    mapping: Mmap,
    _file: File,
}

impl MappedFile {
    /// Maps the complete file as read-only bytes.
    ///
    /// # Errors
    ///
    /// Returns the operating-system error when the file cannot be mapped.
    pub fn new(file: &File) -> io::Result<Self> {
        let retained = file.try_clone()?;
        // SAFETY: the mapping is read-only, the source descriptor is retained
        // for at least as long as the mapping, and this wrapper exposes no
        // operation that can modify or resize the source file.
        let mapping = unsafe { MmapOptions::new().map(file)? };
        Ok(Self {
            mapping,
            _file: retained,
        })
    }

    /// The mapped bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.mapping
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

#[cfg(test)]
#[path = "mapped_tests.rs"]
mod tests;
