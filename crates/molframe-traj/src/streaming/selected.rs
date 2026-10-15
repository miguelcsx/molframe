//! Frame-slice and atom-subset views over any pull reader, and dense chunk filling.
//!
//! [`SelectedReader`] turns "every `step`-th frame in `start..stop`, only these
//! atoms" into an ordinary [`TrajectoryReader`], so the streaming analyses and
//! the chunked Python iterator share one definition of a selection. Sources
//! with random access are sought; forward-only sources decode and discard the
//! frames in between with one reusable scratch buffer, so memory stays
//! independent of trajectory length either way.

use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};
use std::mem::size_of;

/// Every `step`-th frame of `start..stop`, in source order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameSlice {
    /// First source frame.
    pub start: usize,
    /// Exclusive source-frame bound; `usize::MAX` means to the end.
    pub stop: usize,
    /// Distance between selected source frames; never zero.
    pub step: usize,
}

impl FrameSlice {
    /// Every frame of the source.
    pub const ALL: Self = Self {
        start: 0,
        stop: usize::MAX,
        step: 1,
    };

    /// Number of selected frames among `total` source frames.
    #[must_use]
    pub const fn count(&self, total: usize) -> usize {
        let stop = if self.stop < total { self.stop } else { total };
        if self.start >= stop || self.step == 0 {
            0
        } else {
            (stop - self.start).div_ceil(self.step)
        }
    }
}

/// A frame slice and optional atom subset applied to a pull reader.
///
/// Output frames carry the *source* frame ordinal in [`Timestep::frame`] and
/// only positions, times and cells; velocities and forces are not forwarded.
#[derive(Debug)]
pub struct SelectedReader<R> {
    inner: R,
    slice: FrameSlice,
    atoms: Option<Vec<usize>>,
    scratch: Timestep,
    /// Source frame the inner reader will decode next.
    cursor: usize,
    /// Selected frames already emitted.
    emitted: usize,
}

impl<R: TrajectoryReader> SelectedReader<R> {
    /// Wraps `inner`, validating the atom indices against its atom count.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::SelectionOutOfRange`] for an atom index past
    /// the frame, or [`TrajectoryError::InvalidSource`] for a zero stride.
    pub fn new(
        inner: R,
        slice: FrameSlice,
        atoms: Option<Vec<usize>>,
    ) -> Result<Self, TrajectoryError> {
        if slice.step == 0 {
            return Err(TrajectoryError::InvalidSource {
                format: "frame stride must be positive",
            });
        }
        let total = inner.n_atoms();
        if let Some(index) = atoms
            .iter()
            .flatten()
            .copied()
            .find(|index| *index >= total)
        {
            return Err(TrajectoryError::SelectionOutOfRange {
                index,
                atoms: total,
            });
        }
        Ok(Self {
            inner,
            slice,
            atoms,
            scratch: Timestep::default(),
            cursor: 0,
            emitted: 0,
        })
    }

    fn target(&self) -> Option<usize> {
        self.emitted
            .checked_mul(self.slice.step)
            .and_then(|offset| self.slice.start.checked_add(offset))
            .filter(|target| *target < self.slice.stop)
    }

    fn position(
        &mut self,
        target: usize,
        bytes: usize,
        bounded: bool,
    ) -> Result<bool, TrajectoryError> {
        if self.inner.random_access() == RandomAccess::None {
            while self.cursor < target {
                if !decode(&mut self.inner, &mut self.scratch, bytes, bounded)? {
                    return Ok(false);
                }
                self.cursor += 1;
            }
        } else if self.cursor != target {
            self.inner.seek(target)?;
            self.cursor = target;
        }
        Ok(true)
    }

    fn read_selected(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
        bounded: bool,
    ) -> Result<bool, TrajectoryError> {
        let Some(target) = self.target() else {
            return Ok(false);
        };
        if self.inner.n_frames().is_some_and(|total| target >= total)
            || !self.position(target, bytes, bounded)?
        {
            return Ok(false);
        }
        let read = match &self.atoms {
            None => decode(&mut self.inner, timestep, bytes, bounded)?,
            Some(_) => decode(&mut self.inner, &mut self.scratch, bytes, bounded)?,
        };
        if !read {
            return Ok(false);
        }
        self.cursor = target + 1;
        if let Some(atoms) = &self.atoms {
            let source = &self.scratch;
            timestep.positions.clear();
            timestep
                .positions
                .extend(atoms.iter().map(|index| source.positions[*index]));
            timestep.time = source.time;
            timestep.dt = source.dt;
            timestep.cell = source.cell;
        }
        timestep.velocities = None;
        timestep.forces = None;
        timestep.data.clear();
        timestep.frame = target;
        self.emitted += 1;
        Ok(true)
    }

    fn workspace(&self, bytes: usize) -> Result<(), TrajectoryError> {
        let Some(atoms) = &self.atoms else {
            return Ok(());
        };
        let required = self
            .inner
            .n_atoms()
            .checked_add(atoms.len())
            .and_then(|count| count.checked_mul(size_of::<[f32; 3]>()))
            .and_then(|value| value.checked_add(2 * size_of::<Timestep>()))
            .ok_or(TrajectoryError::MemoryLimit {
                required: usize::MAX,
                limit: bytes,
            })?;
        if required > bytes {
            return Err(TrajectoryError::MemoryLimit {
                required,
                limit: bytes,
            });
        }
        Ok(())
    }
}

impl<R: TrajectoryReader> TrajectoryReader for SelectedReader<R> {
    fn format(&self) -> &'static str {
        self.inner.format()
    }

    fn n_atoms(&self) -> usize {
        self.atoms
            .as_ref()
            .map_or_else(|| self.inner.n_atoms(), Vec::len)
    }

    fn n_frames(&self) -> Option<usize> {
        self.inner.n_frames().map(|total| self.slice.count(total))
    }

    fn units(&self) -> Units {
        self.inner.units()
    }

    fn random_access(&self) -> RandomAccess {
        self.inner.random_access()
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.read_selected(timestep, usize::MAX, false)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        self.workspace(bytes)?;
        self.read_selected(timestep, bytes, true)
    }

    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        if self.inner.random_access() == RandomAccess::None {
            let target = frame
                .checked_mul(self.slice.step)
                .and_then(|offset| self.slice.start.checked_add(offset))
                .ok_or(TrajectoryError::RandomAccessUnavailable)?;
            if target < self.cursor {
                return Err(TrajectoryError::RandomAccessUnavailable);
            }
        }
        self.emitted = frame;
        Ok(())
    }

    fn share_index(&mut self) -> Option<crate::FrameIndex> {
        self.inner.share_index()
    }

    fn adopt_index(&mut self, index: &crate::FrameIndex) -> bool {
        self.inner.adopt_index(index)
    }
}

/// A rectangular block of frames filled directly into single buffers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DenseChunk {
    /// Atoms in every frame.
    pub atoms: usize,
    /// Source ordinal of every frame.
    pub indices: Vec<u64>,
    /// `frames * atoms * 3` coordinates, frame then atom then axis.
    pub coordinates: Vec<f32>,
    /// Time of each frame in picoseconds; `NaN` where unrecorded.
    pub times: Vec<f64>,
    /// `frames * 6` cell rows (lengths then angles), `NaN` where a frame has
    /// none; `None` when no frame in the chunk has a cell.
    pub cells: Option<Vec<f64>>,
}

impl DenseChunk {
    /// Appends one frame in place.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::AtomCountMismatch`] when the frame is not
    /// `atoms` wide, or [`TrajectoryError::IdentityOverflow`] for an ordinal
    /// outside 64 bits.
    pub fn push(&mut self, frame: &Timestep) -> Result<(), TrajectoryError> {
        if frame.positions.len() != self.atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.atoms,
                found: frame.positions.len(),
            });
        }
        let index = self.indices.len();
        self.indices
            .push(u64::try_from(frame.frame).map_err(|_| TrajectoryError::IdentityOverflow)?);
        self.coordinates.extend(frame.positions.iter().flatten());
        self.times.push(match frame.time {
            Some(time) => time,
            None => f64::NAN,
        });
        match (frame.cell, &mut self.cells) {
            (Some(cell), Some(rows)) => push_cell(rows, &cell),
            (Some(cell), None) => {
                let mut rows = vec![f64::NAN; index * 6];
                push_cell(&mut rows, &cell);
                self.cells = Some(rows);
            }
            (None, Some(rows)) => rows.extend([f64::NAN; 6]),
            (None, None) => {}
        }
        Ok(())
    }

    /// Number of frames in the chunk.
    #[must_use]
    pub fn frames(&self) -> usize {
        self.indices.len()
    }
}

fn decode<R: TrajectoryReader>(
    reader: &mut R,
    destination: &mut Timestep,
    bytes: usize,
    bounded: bool,
) -> Result<bool, TrajectoryError> {
    if bounded {
        reader.read_next_bounded(destination, bytes)
    } else {
        reader.read_next(destination)
    }
}

fn reserve<T>(vector: &mut Vec<T>, count: usize) -> Result<(), TrajectoryError> {
    vector
        .try_reserve_exact(count)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: count.saturating_mul(size_of::<T>()),
            limit: 0,
        })
}

/// Reads up to `max_frames` frames into freshly allocated, exactly sized buffers.
///
/// Each frame is decoded into `frame` (reused across calls) and appended in
/// place, so peak memory is the output plus one frame, never two copies. The
/// result holds fewer than `max_frames` frames only at end of stream.
///
/// # Errors
///
/// Returns reader errors, or [`TrajectoryError::MemoryLimit`] when the output
/// cannot be allocated.
pub fn read_dense_chunk<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    frame: &mut Timestep,
    max_frames: usize,
) -> Result<DenseChunk, TrajectoryError> {
    let atoms = reader.n_atoms();
    // An open-ended source with an unbounded request grows by amortised doubling
    // instead of reserving a count nobody knows.
    let expected = match reader.n_frames() {
        Some(total) => total.min(max_frames),
        None if max_frames == usize::MAX => 0,
        None => max_frames,
    };
    let scalars = expected
        .checked_mul(atoms)
        .and_then(|count| count.checked_mul(3))
        .ok_or(TrajectoryError::MemoryLimit {
            required: usize::MAX,
            limit: 0,
        })?;
    let mut chunk = DenseChunk {
        atoms,
        ..DenseChunk::default()
    };
    reserve(&mut chunk.coordinates, scalars)?;
    reserve(&mut chunk.times, expected)?;
    reserve(&mut chunk.indices, expected)?;
    while chunk.indices.len() < max_frames && reader.read_next(frame)? {
        chunk.push(frame)?;
    }
    Ok(chunk)
}

fn push_cell(rows: &mut Vec<f64>, cell: &molframe_core::structure::UnitCell) {
    rows.extend(cell.lengths);
    rows.extend(cell.angles);
}

/// Counts the frames a reader yields, decoding each into `frame` and keeping none.
///
/// # Errors
///
/// Returns reader errors.
pub fn count_frames<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    frame: &mut Timestep,
) -> Result<usize, TrajectoryError> {
    let mut count = 0_usize;
    while reader.read_next(frame)? {
        count += 1;
    }
    Ok(count)
}

#[cfg(test)]
#[path = "selected_tests.rs"]
mod tests;
