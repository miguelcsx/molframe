//! Validated borrowed views over fixed-width trajectory coordinates.

use crate::Timestep;

/// Why a borrowed frame view could not be constructed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum FrameViewError {
    /// Frame dimensions cannot be represented as one contiguous coordinate span.
    #[error("frame dimensions overflow addressable coordinates")]
    DimensionOverflow,
    /// The flattened coordinate span is not exactly `frames * atoms` long.
    #[error("coordinate count does not match frame dimensions")]
    DimensionMismatch,
}

/// A validated borrowed view over C-contiguous `(frames, atoms, 3)` coordinates.
///
/// Constructing this view does not copy atom coordinates. Random frame access is
/// constant time and all consumers observe the same backing coordinate span.
#[derive(Clone, Copy, Debug)]
pub struct FrameView<'a> {
    positions: &'a [[f32; 3]],
    frame_count: usize,
    atom_count: usize,
}

impl<'a> FrameView<'a> {
    /// Validates dimensions for one flattened coordinate span.
    ///
    /// # Errors
    ///
    /// Returns [`FrameViewError`] when the dimensions overflow or do not match
    /// the supplied coordinate count.
    pub fn new(
        positions: &'a [[f32; 3]],
        frame_count: usize,
        atom_count: usize,
    ) -> Result<Self, FrameViewError> {
        let expected = frame_count
            .checked_mul(atom_count)
            .ok_or(FrameViewError::DimensionOverflow)?;
        if positions.len() != expected {
            return Err(FrameViewError::DimensionMismatch);
        }
        Ok(Self {
            positions,
            frame_count,
            atom_count,
        })
    }

    /// Number of frames in the view.
    #[must_use]
    pub const fn frame_count(self) -> usize {
        self.frame_count
    }

    /// Number of atoms in every frame.
    #[must_use]
    pub const fn atom_count(self) -> usize {
        self.atom_count
    }

    /// Whether the view has no frames.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.frame_count == 0
    }

    /// Borrows one frame in constant time.
    #[must_use]
    pub fn frame(self, index: usize) -> Option<&'a [[f32; 3]]> {
        if index >= self.frame_count {
            return None;
        }
        let start = index.checked_mul(self.atom_count)?;
        let end = start.checked_add(self.atom_count)?;
        self.positions.get(start..end)
    }
}

/// Internal common access contract for owned and borrowed frame sources.
pub(crate) trait FrameSource {
    /// Number of frames available to a kernel.
    fn frame_count(&self) -> usize;
    /// Number of atoms in the first frame, or zero for an empty source.
    fn atom_count(&self) -> usize;
    /// One frame in coordinate order.
    fn frame(&self, index: usize) -> Option<&[[f32; 3]]>;
}

impl FrameSource for FrameView<'_> {
    fn frame_count(&self) -> usize {
        FrameView::frame_count(*self)
    }

    fn atom_count(&self) -> usize {
        FrameView::atom_count(*self)
    }

    fn frame(&self, index: usize) -> Option<&[[f32; 3]]> {
        (*self).frame(index)
    }
}

impl FrameSource for [Vec<[f32; 3]>] {
    fn frame_count(&self) -> usize {
        self.len()
    }

    fn atom_count(&self) -> usize {
        match self.first() {
            Some(frame) => frame.len(),
            None => 0,
        }
    }

    fn frame(&self, index: usize) -> Option<&[[f32; 3]]> {
        self.get(index).map(Vec::as_slice)
    }
}

impl FrameSource for [Timestep] {
    fn frame_count(&self) -> usize {
        self.len()
    }

    fn atom_count(&self) -> usize {
        match self.first() {
            Some(frame) => frame.positions.len(),
            None => 0,
        }
    }

    fn frame(&self, index: usize) -> Option<&[[f32; 3]]> {
        self.get(index).map(|frame| frame.positions.as_slice())
    }
}

#[cfg(test)]
#[path = "frame_view_tests.rs"]
mod tests;
