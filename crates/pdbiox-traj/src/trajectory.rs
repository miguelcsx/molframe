//! The in-memory trajectory: frames with random access and slicing.
//!
//! Holding every frame in memory is the simplest source — an array of coordinate
//! sets — and it gives constant-time random access to any frame and a cheap way
//! to take a subset by index. Streaming formats that cannot afford to hold every
//! frame are a separate concern; this is the model those readers, and callers
//! with coordinates already in hand, hand their frames to.

/// One frame: the atom positions at a single point in a trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    /// The atom positions, one per atom, in atom order.
    pub positions: Vec<[f32; 3]>,
}

/// A sequence of frames over a fixed set of atoms.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trajectory {
    frames: Vec<Frame>,
}

impl Trajectory {
    /// Builds a trajectory from frames already in memory.
    #[must_use]
    pub fn from_frames(frames: Vec<Frame>) -> Self {
        Self { frames }
    }

    /// The number of frames.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether the trajectory has no frames.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// The frame at `index`, or `None` when it is out of range.
    #[must_use]
    pub fn frame(&self, index: usize) -> Option<&Frame> {
        self.frames.get(index)
    }

    /// All frames in order.
    #[must_use]
    pub fn frames(&self) -> &[Frame] {
        &self.frames
    }

    /// A new trajectory of the frames at the given indices, in the given order.
    ///
    /// Out-of-range indices are skipped, so slicing never fails; an index may be
    /// repeated to duplicate a frame.
    #[must_use]
    pub fn slice(&self, indices: &[usize]) -> Trajectory {
        let frames = indices
            .iter()
            .filter_map(|&index| self.frames.get(index).cloned())
            .collect();
        Trajectory { frames }
    }
}

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
