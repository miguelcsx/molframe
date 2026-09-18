//! Text scalar grids: Gaussian cube and `OpenDX`.
//!
//! Both formats describe the same thing an MRC map does — a scalar sampled on
//! a regular lattice — so both land on [`DensityMap`] and everything already
//! built on that (statistics, histograms, Cartesian sampling) works unchanged.
//!
//! Two things differ from MRC and are handled here rather than left to the
//! caller. Both formats store the *last* axis fastest, where a map stores the
//! first, so parsed values are placed directly at their canonical index; and a
//! cube's geometry is conventionally in Bohr, so its axes are converted to
//! angstroms. Reading is `O(values)` with one value-sized allocation rather
//! than retaining both file-order and canonical volumes.

#[path = "grid/cube.rs"]
mod cube;
#[path = "grid/dx.rs"]
mod dx;

pub use cube::{CubeAtom, CubeGrid, read_cube};
pub use dx::read_dx;

use molframe_core::structure::UnitCell;

use crate::mrc::DensityMap;

/// One bohr in angstroms.
///
/// Cube geometry is quantum-chemistry native and therefore atomic units; the
/// rest of the library is angstroms, and converting at the boundary keeps every
/// downstream distance comparable.
pub(crate) const BOHR_ANGSTROMS: f64 = 0.529_177_210_903;

/// Largest grid a text reader will reserve for, in values.
///
/// A malformed header can claim any dimensions it likes. Refusing an
/// implausible product up front means a corrupt file costs a parse error rather
/// than the host's memory.
pub(crate) const MAX_GRID_VALUES: usize = 1 << 30;

/// Invalid or unsupported text grid data.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum GridError {
    /// The file ended before its declared content.
    #[error("truncated {format} grid")]
    Truncated {
        /// Format name for the diagnostic.
        format: &'static str,
    },
    /// A header line was missing, malformed, or out of order.
    #[error("invalid {format} header: {reason}")]
    InvalidHeader {
        /// Format name for the diagnostic.
        format: &'static str,
        /// What was expected at that point.
        reason: &'static str,
    },
    /// A number could not be read where one was required.
    #[error("invalid {format} number at value {index}")]
    InvalidNumber {
        /// Format name for the diagnostic.
        format: &'static str,
        /// Position in the value stream.
        index: usize,
    },
    /// Declared dimensions overflow or exceed the reader's ceiling.
    #[error("{format} grid dimensions are out of range")]
    SizeOverflow {
        /// Format name for the diagnostic.
        format: &'static str,
    },
    /// The host could not reserve the validated value buffer.
    #[error("{format} grid allocation exceeds host resources")]
    ResourceLimit {
        /// Format name for the diagnostic.
        format: &'static str,
    },
    /// The text is not UTF-8.
    #[error("{format} grid is not valid UTF-8")]
    NotUtf8 {
        /// Format name for the diagnostic.
        format: &'static str,
    },
}

/// The three lattice vectors spanning one voxel, in angstroms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Axes {
    pub(crate) origin: [f64; 3],
    pub(crate) steps: [[f64; 3]; 3],
    pub(crate) counts: [usize; 3],
}

impl Axes {
    /// The unit cell the whole grid spans.
    ///
    /// Cube and DX both describe geometry as per-voxel step vectors, which is
    /// the same information a cell carries once each is multiplied by its count
    /// — so a grid whose axes are skewed survives the conversion instead of
    /// being silently forced orthogonal.
    pub(crate) fn cell(&self) -> UnitCell {
        let spans = self.spans();
        let lengths = spans.map(norm);
        let angles = [
            angle_between(spans[1], spans[2]),
            angle_between(spans[0], spans[2]),
            angle_between(spans[0], spans[1]),
        ];
        UnitCell { lengths, angles }
    }

    fn spans(&self) -> [[f64; 3]; 3] {
        let mut spans = [[0.0; 3]; 3];
        for (axis, (span, step)) in spans.iter_mut().zip(self.steps).enumerate() {
            let count = usize_to_f64(self.counts[axis]);
            for (value, component) in span.iter_mut().zip(step) {
                *value = component * count;
            }
        }
        spans
    }
}

/// Maps a scalar from last-axis-fastest file order to first-axis-fastest map
/// order while preserving adjacent fields within each voxel.
pub(crate) fn canonical_value_index(source: usize, counts: [usize; 3], fields: usize) -> usize {
    let field = source % fields;
    let voxel = source / fields;
    let z = voxel % counts[2];
    let xy = voxel / counts[2];
    let y = xy % counts[1];
    let x = xy / counts[1];
    (((z * counts[1] + y) * counts[0] + x) * fields) + field
}

/// Assembles the canonical map once its geometry and values are known.
pub(crate) fn density_map(axes: &Axes, values: Vec<f32>) -> DensityMap {
    DensityMap {
        dimensions: axes.counts,
        starts: [0; 3],
        sampling: axes.counts,
        cell: axes.cell(),
        origin: axes.origin,
        space_group: 0,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values,
    }
}

/// Multiplies the three counts, refusing anything past the reader's ceiling.
pub(crate) fn value_count(counts: [usize; 3], format: &'static str) -> Result<usize, GridError> {
    let product = counts
        .iter()
        .try_fold(1_usize, |product, count| product.checked_mul(*count));
    match product {
        Some(product) if product <= MAX_GRID_VALUES && product > 0 => Ok(product),
        _ => Err(GridError::SizeOverflow { format }),
    }
}

/// Reserves exactly the validated buffer, reporting a refusal rather than
/// aborting when the host cannot supply it.
pub(crate) fn reserve_values(count: usize, format: &'static str) -> Result<Vec<f32>, GridError> {
    let mut values = Vec::new();
    match values.try_reserve_exact(count) {
        Ok(()) => Ok(values),
        Err(_) => Err(GridError::ResourceLimit { format }),
    }
}

/// Allocates one complete output volume so parsers can scatter directly into
/// canonical order without retaining a second file-order volume.
pub(crate) fn zeroed_values(count: usize, format: &'static str) -> Result<Vec<f32>, GridError> {
    let mut values = reserve_values(count, format)?;
    values.resize(count, 0.0);
    Ok(values)
}

fn norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|value| value * value).sum::<f64>().sqrt()
}

fn angle_between(left: [f64; 3], right: [f64; 3]) -> f64 {
    let denominator = norm(left) * norm(right);
    if denominator == 0.0 {
        return 90.0;
    }
    let cosine = left
        .iter()
        .zip(right)
        .map(|(first, second)| first * second)
        .sum::<f64>()
        / denominator;
    cosine.clamp(-1.0, 1.0).acos().to_degrees()
}

fn usize_to_f64(value: usize) -> f64 {
    crate::numeric::usize_to_f64(value)
}

#[cfg(test)]
#[path = "grid_tests.rs"]
mod tests;
