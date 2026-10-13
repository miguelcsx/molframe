//! Open pockets found by scanning a grid, after LIGSITE (Hendlich et al. 1997): a grid
//! point outside the molecule is *buried* when, along most of seven scan directions
//! (three axes, four body diagonals), atoms lie on both sides of it within reach.
//! Connected buried points form a pocket; pockets are returned largest first.
//!
//! This finds geometric concavities, nothing more. It does not know chemistry,
//! druggability or whether a ligand binds there. [`cavities`](crate::cavities) reports
//! only sealed cavities, which are not binding pockets, hence this. A buried point must
//! also be reachable by a probe from outside, which leaves out packing voids and sealed
//! cavities.
//!
//! The painting and the scan run over slabs of the grid in parallel and give the cells a
//! serial pass would; the flood and the labelling are linear in the number of cells.
//! Cancellation is polled between and inside the passes.

use crate::numeric::{f32_to_usize, usize_to_f32};
use molframe_core::{CancellationToken, Cancelled};
use thiserror::Error;

#[path = "pocket_grid.rs"]
mod grid;
use grid::{Grid, buried, components, paint, reachable};

#[path = "pocket_points.rs"]
mod points;
pub use points::PointProximity;

/// The parameters of the pocket scan. [`Default`] is the published LIGSITE setting at
/// one ångström spacing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PocketScanOptions {
    /// Grid spacing, ångström.
    pub spacing: f32,
    /// How far a scan looks for atoms on each side, ångström.
    pub reach: f32,
    /// A point is buried when this many of the seven directions find atoms on both sides.
    pub min_directions: usize,
    /// Pockets with fewer grid points than this are noise.
    pub min_points: usize,
    /// A pocket point must be somewhere a probe of this radius can reach from outside,
    /// ångström.
    pub probe: f32,
    /// Empty space kept around the atoms' box, ångström.
    pub margin: f32,
    /// The most grid cells the scan accepts.
    pub max_cells: usize,
    /// Fewer atoms than this have no pocket to find.
    pub min_atoms: usize,
}

impl Default for PocketScanOptions {
    fn default() -> Self {
        Self {
            spacing: 1.0,
            reach: 12.0,
            min_directions: 5,
            min_points: 24,
            probe: 1.2,
            margin: 4.0,
            max_cells: 40_000_000,
            min_atoms: 10,
        }
    }
}

/// One pocket: the grid points that form it.
#[derive(Clone, Debug, PartialEq)]
pub struct GridPocket {
    /// The grid points, ångström, in scan order.
    pub points: Vec<[f32; 3]>,
    /// The mean of the points, ångström.
    pub centre: [f32; 3],
    /// The volume of the points, cubic ångström.
    pub volume: f32,
}

/// Why a scan produced no answer; "no pocket" would be a claim the scan did not make.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum PocketScanError {
    /// Positions and radii differ in length.
    #[error("{positions} positions but {radii} radii")]
    LengthMismatch {
        /// Positions given.
        positions: usize,
        /// Radii given.
        radii: usize,
    },
    /// The spacing, reach or probe is not a finite positive (or zero probe) number.
    #[error("the pocket scan options must be finite, with a positive spacing")]
    InvalidOptions,
    /// Fewer atoms than the scan needs to hold a pocket.
    #[error("only {atoms} atoms: too few to hold a pocket")]
    TooFewAtoms {
        /// How many there were.
        atoms: usize,
    },
    /// The grid that would hold the molecule is larger than the scan accepts.
    #[error("the grid would have {cells} cells, more than the {limit} the scan takes on")]
    TooLarge {
        /// Cells the grid would have.
        cells: usize,
        /// The most the scan accepts.
        limit: usize,
    },
    /// The scan was cancelled before it finished.
    #[error("the pocket scan was cancelled")]
    Cancelled,
}

impl From<Cancelled> for PocketScanError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// The pockets of atoms at `positions` with van der Waals `radii`, largest first, at
/// most `limit`.
///
/// # Errors
///
/// [`PocketScanError`] for mismatched or invalid input, a molecule too small to hold a
/// pocket or too large to scan, or when `cancel` was triggered.
pub fn scan_pockets(
    positions: &[[f32; 3]],
    radii: &[f32],
    limit: usize,
    options: &PocketScanOptions,
    cancel: &CancellationToken,
) -> Result<Vec<GridPocket>, PocketScanError> {
    if positions.len() != radii.len() {
        return Err(PocketScanError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }
    let finite = [
        options.spacing,
        options.reach,
        options.probe,
        options.margin,
    ]
    .iter()
    .all(|v| v.is_finite());
    if !finite || options.spacing <= 0.0 {
        return Err(PocketScanError::InvalidOptions);
    }
    if positions.len() < options.min_atoms {
        return Err(PocketScanError::TooFewAtoms {
            atoms: positions.len(),
        });
    }
    let (mut lo, mut hi) = ([f32::MAX; 3], [f32::MIN; 3]);
    for p in positions {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let spacing = options.spacing;
    let dims = [0, 1, 2].map(|k| {
        f32_to_usize(((hi[k] - lo[k] + 2.0 * options.margin) / spacing).ceil()).saturating_add(1)
    });
    let cells = dims[0].saturating_mul(dims[1]).saturating_mul(dims[2]);
    if cells > options.max_cells {
        // a molecule this large is not a grid to scan blindly
        return Err(PocketScanError::TooLarge {
            cells,
            limit: options.max_cells,
        });
    }
    let grid = Grid {
        origin: lo.map(|c| c - options.margin),
        dims,
        spacing,
    };
    // Atom cells: inside an atom's sphere. Probe-blocked cells: inside a sphere grown by
    // the probe, where a probe's centre cannot be.
    let solid = paint(&grid, positions, radii, 0.0, cancel)?;
    let blocked = paint(&grid, positions, radii, options.probe, cancel)?;
    let reachable = reachable(&grid, &blocked, cancel)?;
    drop(blocked);
    let buried = buried(&grid, &solid, &reachable, options, cancel)?;
    let mut groups = components(&grid, &buried, options.min_points, cancel)?;
    groups.sort_by_key(|g| std::cmp::Reverse(g.len()));
    groups.truncate(limit);
    Ok(groups
        .into_iter()
        .map(|members| {
            let points: Vec<[f32; 3]> = members
                .iter()
                .map(|m| grid.point(m[0], m[1], m[2]))
                .collect();
            let count = usize_to_f32(points.len());
            let centre = [0, 1, 2].map(|k| points.iter().map(|p| p[k]).sum::<f32>() / count);
            GridPocket {
                volume: count * spacing.powi(3),
                points,
                centre,
            }
        })
        .collect())
}

#[cfg(test)]
#[path = "pocket_tests.rs"]
mod tests;
