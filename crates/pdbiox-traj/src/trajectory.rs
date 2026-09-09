//! Contiguous in-memory trajectories with constant-time borrowed frame views.
//!
//! Coordinates use one flat `(frames, atoms, 3)` allocation. Random frame
//! access is `O(1)`, slicing is `O(selected frames * atoms)`, and no frame owns
//! a second allocation. Streaming readers remain the bounded source for data
//! that should not be materialised.

use crate::{FrameView, FrameViewError};

/// One owned frame accepted at an API boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    /// Atom positions in topology order.
    pub positions: Vec<[f32; 3]>,
}

/// One borrowed frame over contiguous trajectory storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameRef<'a> {
    /// Atom positions in topology order.
    pub positions: &'a [[f32; 3]],
}

/// Failure to construct one fixed-width trajectory.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum TrajectoryBuildError {
    /// Coordinate dimensions exceed the addressable process range.
    #[error("trajectory coordinate dimensions overflow addressable memory")]
    DimensionOverflow,
    /// A frame disagrees with the atom count established by the first frame.
    #[error("trajectory frame {frame} has {found} atoms; expected {expected}")]
    AtomCountMismatch {
        /// Zero-based frame index.
        frame: usize,
        /// Fixed atom count established by the first frame.
        expected: usize,
        /// Atom count in the rejected frame.
        found: usize,
    },
}

/// A fixed-width trajectory stored in one contiguous coordinate allocation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trajectory {
    positions: Box<[[f32; 3]]>,
    frames: usize,
    atoms: usize,
}

impl Trajectory {
    /// Flattens fixed-width frames into one coordinate allocation.
    ///
    /// # Errors
    ///
    /// Returns a typed error for ragged frames or dimension overflow.
    pub fn from_frames(frames: Vec<Frame>) -> Result<Self, TrajectoryBuildError> {
        let frame_count = frames.len();
        let atom_count = match frames.first() {
            Some(frame) => frame.positions.len(),
            None => 0,
        };
        let total = frame_count
            .checked_mul(atom_count)
            .ok_or(TrajectoryBuildError::DimensionOverflow)?;
        for (index, frame) in frames.iter().enumerate() {
            if frame.positions.len() != atom_count {
                return Err(TrajectoryBuildError::AtomCountMismatch {
                    frame: index,
                    expected: atom_count,
                    found: frame.positions.len(),
                });
            }
        }
        let mut positions = Vec::with_capacity(total);
        for frame in frames {
            positions.extend(frame.positions);
        }
        Ok(Self {
            positions: positions.into_boxed_slice(),
            frames: frame_count,
            atoms: atom_count,
        })
    }

    /// Number of frames.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.frames
    }

    /// Number of atoms in every frame.
    #[must_use]
    pub const fn atom_count(&self) -> usize {
        self.atoms
    }

    /// Whether the trajectory has no frames.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.frames == 0
    }

    /// Borrows one frame in constant time.
    #[must_use]
    pub fn frame(&self, index: usize) -> Option<FrameRef<'_>> {
        if index >= self.frames {
            return None;
        }
        let start = index.checked_mul(self.atoms)?;
        let end = start.checked_add(self.atoms)?;
        self.positions
            .get(start..end)
            .map(|positions| FrameRef { positions })
    }

    /// Borrows the complete flattened coordinate span.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        &self.positions
    }

    /// Borrows this trajectory through the shared fixed-width frame contract.
    ///
    /// # Errors
    ///
    /// Returns a dimension error if internal invariants are violated.
    pub fn as_view(&self) -> Result<FrameView<'_>, FrameViewError> {
        FrameView::new(&self.positions, self.frames, self.atoms)
    }

    /// Copies selected frames into one new contiguous trajectory.
    ///
    /// Out-of-range indices are skipped and repeated indices duplicate frames.
    /// # Errors
    ///
    /// Returns a typed error when the selected coordinate dimensions overflow.
    pub fn slice(&self, indices: &[usize]) -> Result<Self, TrajectoryBuildError> {
        let retained = indices.iter().filter(|&&index| index < self.frames).count();
        let capacity = retained
            .checked_mul(self.atoms)
            .ok_or(TrajectoryBuildError::DimensionOverflow)?;
        let mut positions = Vec::with_capacity(capacity);
        for &index in indices {
            if let Some(frame) = self.frame(index) {
                positions.extend_from_slice(frame.positions);
            }
        }
        Ok(Self {
            positions: positions.into_boxed_slice(),
            frames: retained,
            atoms: self.atoms,
        })
    }
}

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
