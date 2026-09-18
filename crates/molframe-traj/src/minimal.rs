//! Atom-count-only topology for coordinate formats without chemical metadata.

use crate::{TrajectoryError, TrajectoryReader};

/// The smallest topology that can safely pair with a trajectory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinimalTopology {
    atom_count: usize,
}

impl MinimalTopology {
    /// Declares the fixed number of atoms represented by every frame.
    #[must_use]
    pub const fn new(atom_count: usize) -> Self {
        Self { atom_count }
    }

    /// Declared atom count.
    #[must_use]
    pub const fn atom_count(self) -> usize {
        self.atom_count
    }

    /// Checks the topology before pairing it with a coordinate reader.
    ///
    /// # Errors
    ///
    /// Returns an atom-count mismatch without consuming a frame.
    pub fn validate_reader<R: TrajectoryReader>(self, reader: &R) -> Result<(), TrajectoryError> {
        let found = reader.n_atoms();
        if found != self.atom_count {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.atom_count,
                found,
            });
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "minimal_tests.rs"]
mod tests;
