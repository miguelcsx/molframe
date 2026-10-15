//! Several homogeneous readers presented as one forward sequence.

use super::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};

/// Several homogeneous readers presented as one forward sequence.
#[derive(Debug)]
pub struct ChainedReader<R> {
    readers: Vec<R>,
    current: usize,
    frame: usize,
    atoms: usize,
}

impl<R: TrajectoryReader> ChainedReader<R> {
    /// Validates atom counts before consuming any source.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::AtomCountMismatch`] if sources cannot share a
    /// topology.
    pub fn new(readers: Vec<R>) -> Result<Self, TrajectoryError> {
        let atoms = readers.first().map_or(0, TrajectoryReader::n_atoms);
        for reader in &readers {
            if reader.n_atoms() != atoms {
                return Err(TrajectoryError::AtomCountMismatch {
                    expected: atoms,
                    found: reader.n_atoms(),
                });
            }
        }
        Ok(Self {
            readers,
            current: 0,
            frame: 0,
            atoms,
        })
    }
}

impl<R: TrajectoryReader> TrajectoryReader for ChainedReader<R> {
    fn format(&self) -> &'static str {
        "chain"
    }

    fn n_atoms(&self) -> usize {
        self.atoms
    }

    fn n_frames(&self) -> Option<usize> {
        self.readers.iter().try_fold(0usize, |total, reader| {
            reader.n_frames().and_then(|count| total.checked_add(count))
        })
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        RandomAccess::None
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.next_from_sources(timestep, TrajectoryReader::read_next)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        self.next_from_sources(timestep, |reader, timestep| {
            reader.read_next_bounded(timestep, bytes)
        })
    }

    fn seek(&mut self, _frame: usize) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::RandomAccessUnavailable)
    }
}

impl<R: TrajectoryReader> ChainedReader<R> {
    /// Reads from the current source, moving on when one is exhausted, and
    /// renumbers the frame across the whole chain.
    fn next_from_sources(
        &mut self,
        timestep: &mut Timestep,
        mut read: impl FnMut(&mut R, &mut Timestep) -> Result<bool, TrajectoryError>,
    ) -> Result<bool, TrajectoryError> {
        while let Some(reader) = self.readers.get_mut(self.current) {
            if read(reader, timestep)? {
                timestep.frame = self.frame;
                self.frame += 1;
                return Ok(true);
            }
            self.current += 1;
        }
        Ok(false)
    }
}
