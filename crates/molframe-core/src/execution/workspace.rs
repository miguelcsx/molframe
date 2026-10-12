//! Reserved private directories for native tools. Files are checked in O(files).

use super::spill::DiskAccount;
use super::{ExecutionContext, SpillError};
use std::{fs, io, path::Path, sync::Arc};

/// Private native-tool files retaining one shared temporary-storage reservation.
///
/// External writers cannot be prevented portably from exceeding a reservation.
/// Hosts must call `check_size` while the tool runs, terminate an offending
/// writer, and check again before consuming its outputs.
#[derive(Debug)]
pub struct TemporaryWorkspace {
    directory: tempfile::TempDir,
    account: Arc<DiskAccount>,
    reserved: u64,
    maximum_files: usize,
}

impl ExecutionContext {
    /// Reserves temporary bytes before creating one private native-tool directory.
    /// Concurrent tools and spill files share this context's disk ceiling.
    ///
    /// # Errors
    /// Returns disabled storage, exhausted reservation or filesystem errors.
    pub fn temporary_workspace(
        &self,
        maximum_bytes: u64,
    ) -> Result<TemporaryWorkspace, SpillError> {
        self.temporary_workspace_with_file_limit(maximum_bytes, 64)
    }

    /// Reserves a private directory with an explicit regular-file count ceiling.
    ///
    /// # Errors
    /// Refuses zero file capacity, disabled storage or an exhausted disk budget.
    pub fn temporary_workspace_with_file_limit(
        &self,
        maximum_bytes: u64,
        maximum_files: usize,
    ) -> Result<TemporaryWorkspace, SpillError> {
        if maximum_files == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "workspace file capacity must be positive",
            )
            .into());
        }
        let Some(account) = &self.disk else {
            return Err(SpillError::Disabled);
        };
        if maximum_bytes == 0 {
            return Err(SpillError::BudgetExceeded {
                requested: 1,
                available: 0,
            });
        }
        account.try_charge(maximum_bytes)?;
        let directory = fs::create_dir_all(&account.root).and_then(|()| {
            tempfile::Builder::new()
                .prefix("molframe-tool-")
                .tempdir_in(&account.root)
        });
        let directory = match directory {
            Ok(directory) => directory,
            Err(error) => {
                account.release(maximum_bytes);
                return Err(error.into());
            }
        };
        Ok(TemporaryWorkspace {
            directory,
            account: Arc::clone(account),
            reserved: maximum_bytes,
            maximum_files,
        })
    }
}

impl TemporaryWorkspace {
    /// Directory owned by this reservation, never a caller's existing directory.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.directory.path()
    }

    /// Checks regular files against the retained byte and file-count reservations.
    /// Subdirectories and symlinks are refused rather than traversed.
    ///
    /// # Errors
    /// Returns a disk-budget or I/O error for oversized or unsupported output.
    pub fn check_size(&self) -> Result<u64, SpillError> {
        let mut bytes = 0u64;
        for (index, entry) in fs::read_dir(self.path())?.enumerate() {
            let entry = entry?;
            if index >= self.maximum_files || !entry.file_type()?.is_file() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "native workspace file capacity or regular-file requirement exceeded",
                )
                .into());
            }
            bytes = bytes.saturating_add(entry.metadata()?.len());
            if bytes > self.reserved {
                return Err(SpillError::BudgetExceeded {
                    requested: bytes,
                    available: self.reserved,
                });
            }
        }
        Ok(bytes)
    }
}

impl Drop for TemporaryWorkspace {
    fn drop(&mut self) {
        match fs::remove_dir_all(self.path()) {
            Ok(()) => self.account.release(self.reserved),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                self.account.release(self.reserved);
            }
            // Failed cleanup keeps its charge: those files may still occupy disk.
            Err(_) => {}
        }
    }
}

#[cfg(test)]
#[path = "workspace_tests.rs"]
mod tests;
