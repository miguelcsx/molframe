//! Borrowed data a plan execution reads from.
//!
//! Every input is a view over memory the caller owns. Nothing here allocates or
//! copies a coordinate, index, scalar or mask array, which is what lets a
//! binding drive a plan straight from a `NumPy` buffer or a C-contiguous block
//! without an intermediate representation.

use molframe_core::structure::Structure;

/// One borrowed coordinate array supplied to a plan execution.
#[derive(Clone, Copy, Debug)]
pub struct CoordinateInput<'a> {
    /// C-contiguous coordinate rows. The facade never materialises this slice.
    pub positions: &'a [[f32; 3]],
}

/// A borrowed flattened view over a C-contiguous `(frames, atoms, 3)` array.
#[derive(Clone, Copy, Debug)]
pub struct FrameInput<'a> {
    /// Flattened coordinate rows, retained by the binding owner.
    pub positions: &'a [[f32; 3]],
    /// Number of trajectory frames represented by `positions`.
    pub frame_count: usize,
    /// Number of atoms in each frame.
    pub atom_count: usize,
}

/// A borrowed C-contiguous atom-index array supplied to a plan execution.
#[derive(Clone, Copy, Debug)]
pub struct IndexInput<'a> {
    /// Zero-based atom indices. The facade never materialises this slice.
    pub indices: &'a [usize],
}

/// A borrowed one-dimensional `f64` input such as per-atom masses.
#[derive(Clone, Copy, Debug)]
pub struct ScalarInput<'a> {
    /// C-contiguous scalar values.
    pub values: &'a [f64],
}

/// A borrowed one-dimensional `f32` input used by atom-aligned analyses.
#[derive(Clone, Copy, Debug)]
pub struct FloatInput<'a> {
    /// C-contiguous floating-point values.
    pub values: &'a [f32],
}

/// A borrowed one-dimensional boolean input used by surface operations.
#[cfg(feature = "surface")]
#[derive(Clone, Copy, Debug)]
pub struct MaskInput<'a> {
    /// C-contiguous boolean values.
    pub values: &'a [bool],
}

/// Inputs shared by a plan execution.
#[derive(Clone, Copy, Debug, Default)]
pub struct PlanInput<'a> {
    /// Structure snapshot used by structure-bound operations.
    pub structure: Option<&'a Structure>,
    /// Borrowed arrays used by array-bound operations.
    pub arrays: &'a [CoordinateInput<'a>],
    /// Borrowed scalar arrays used by weighted geometry operations.
    pub scalars: &'a [ScalarInput<'a>],
    /// Borrowed floating-point arrays used by atom-aligned analyses.
    pub floats: &'a [FloatInput<'a>],
    /// Borrowed boolean arrays used by surface operations.
    #[cfg(feature = "surface")]
    pub masks: &'a [MaskInput<'a>],
    /// Borrowed frame arrays used by ensemble geometry operations.
    pub frames: &'a [FrameInput<'a>],
    /// Borrowed atom-index arrays used by trajectory analyses.
    pub indices: &'a [IndexInput<'a>],
}
