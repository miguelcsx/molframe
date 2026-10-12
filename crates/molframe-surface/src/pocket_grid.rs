//! The grid passes of the pocket scan: painting atoms into cells, flooding from outside,
//! scanning for buriedness and labelling the connected cells. The first three run over
//! slabs of the first grid axis in parallel and give the cells a serial pass would.

use super::PocketScanOptions;
use crate::numeric::{f32_to_i32, i32_to_usize, usize_to_f32, usize_to_i32};
use molframe_core::{CancellationToken, Cancelled};
use rayon::prelude::*;

/// The scan directions: three axes and four body diagonals.
const DIRECTIONS: [[i32; 3]; 7] = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 1, 1],
    [1, 1, -1],
    [1, -1, 1],
    [-1, 1, 1],
];

/// The face neighbours a flood steps through.
const FACE_STEPS: [[i32; 3]; 6] = [
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
];

pub(super) struct Grid {
    pub(super) origin: [f32; 3],
    pub(super) dims: [usize; 3],
    pub(super) spacing: f32,
}

impl Grid {
    pub(super) fn index(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.dims[1] + j) * self.dims[2] + k
    }

    pub(super) fn point(&self, i: usize, j: usize, k: usize) -> [f32; 3] {
        [
            self.origin[0] + usize_to_f32(i) * self.spacing,
            self.origin[1] + usize_to_f32(j) * self.spacing,
            self.origin[2] + usize_to_f32(k) * self.spacing,
        ]
    }

    pub(super) fn cells(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Cells in one slab of the first axis.
    fn slab(&self) -> usize {
        (self.dims[1] * self.dims[2]).max(1)
    }

    /// The cell `step` away from `from`, when it is inside the grid.
    fn neighbour(&self, from: [usize; 3], step: [i32; 3]) -> Option<[usize; 3]> {
        let mut out = [0usize; 3];
        for axis in 0..3 {
            let at = usize_to_i32(from[axis]) + step[axis];
            let at = usize::try_from(at).ok()?;
            if at >= self.dims[axis] {
                return None;
            }
            out[axis] = at;
        }
        Some(out)
    }

    fn in_cell(&self, cell: [usize; 3]) -> usize {
        self.index(cell[0], cell[1], cell[2])
    }
}

/// An atom's footprint on the grid: its centre cell, how many cells it reaches, and its
/// grown radius.
struct Footprint {
    centre: [i32; 3],
    reach: i32,
    radius: f32,
}

/// Cells inside a sphere of each atom grown by `pad`.
pub(super) fn paint(
    grid: &Grid,
    positions: &[[f32; 3]],
    radii: &[f32],
    pad: f32,
    cancel: &CancellationToken,
) -> Result<Vec<bool>, Cancelled> {
    let spacing = grid.spacing;
    let footprints: Vec<Footprint> = positions
        .iter()
        .zip(radii)
        .map(|(p, r)| {
            let radius = r + pad;
            Footprint {
                centre: [0, 1, 2].map(|k| f32_to_i32(((p[k] - grid.origin[k]) / spacing).round())),
                reach: f32_to_i32((radius / spacing).ceil()),
                radius,
            }
        })
        .collect();
    // Which atoms reach each slab, so a slab only sees its own.
    let mut by_slab: Vec<Vec<usize>> = vec![Vec::new(); grid.dims[0]];
    let last = usize_to_i32(grid.dims[0]) - 1;
    for (n, f) in footprints.iter().enumerate() {
        let from = (f.centre[0] - f.reach).max(0);
        let to = (f.centre[0] + f.reach).min(last);
        for i in from..=to {
            if let Some(slab) = by_slab.get_mut(i32_to_usize(i)) {
                slab.push(n);
            }
        }
    }
    let mut mask = vec![false; grid.cells()];
    mask.par_chunks_mut(grid.slab())
        .enumerate()
        .for_each(|(i, slab)| {
            if cancel.is_cancelled() {
                return;
            }
            for &n in &by_slab[i] {
                let f = &footprints[n];
                let p = positions[n];
                for j in (f.centre[1] - f.reach)..=(f.centre[1] + f.reach) {
                    let Ok(ju) = usize::try_from(j) else { continue };
                    if ju >= grid.dims[1] {
                        continue;
                    }
                    for k in (f.centre[2] - f.reach)..=(f.centre[2] + f.reach) {
                        let Ok(ku) = usize::try_from(k) else { continue };
                        if ku >= grid.dims[2] {
                            continue;
                        }
                        let q = grid.point(i, ju, ku);
                        let d2 =
                            (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2);
                        if d2 <= f.radius * f.radius {
                            slab[ju * grid.dims[2] + ku] = true;
                        }
                    }
                }
            }
        });
    cancel.check()?;
    Ok(mask)
}

/// The cells a probe can reach from outside: flooded in from the grid's faces through
/// cells that are not blocked.
pub(super) fn reachable(
    grid: &Grid,
    blocked: &[bool],
    cancel: &CancellationToken,
) -> Result<Vec<bool>, Cancelled> {
    let mut seen = vec![false; grid.cells()];
    let mut stack = Vec::new();
    for i in 0..grid.dims[0] {
        for j in 0..grid.dims[1] {
            for k in 0..grid.dims[2] {
                let edge = i == 0
                    || j == 0
                    || k == 0
                    || i + 1 == grid.dims[0]
                    || j + 1 == grid.dims[1]
                    || k + 1 == grid.dims[2];
                let at = grid.index(i, j, k);
                if edge && !blocked[at] {
                    seen[at] = true;
                    stack.push([i, j, k]);
                }
            }
        }
    }
    let mut steps = 0usize;
    while let Some(cell) = stack.pop() {
        steps += 1;
        if steps.is_multiple_of(1 << 20) {
            cancel.check()?;
        }
        for step in FACE_STEPS {
            let Some(next) = grid.neighbour(cell, step) else {
                continue;
            };
            let n = grid.in_cell(next);
            if !blocked[n] && !seen[n] {
                seen[n] = true;
                stack.push(next);
            }
        }
    }
    Ok(seen)
}

/// Cells outside the protein, reachable from outside, with protein on both sides along at
/// least `options.min_directions` of the seven scan directions.
pub(super) fn buried(
    grid: &Grid,
    solid: &[bool],
    reachable: &[bool],
    options: &PocketScanOptions,
    cancel: &CancellationToken,
) -> Result<Vec<bool>, Cancelled> {
    let steps = f32_to_i32(options.reach / grid.spacing);
    let hits = |cell: [usize; 3], d: [i32; 3], sign: i32| -> bool {
        let step = if d.iter().filter(|c| **c != 0).count() == 1 {
            1.0
        } else {
            3f32.sqrt()
        };
        let reach = f32_to_i32(options.reach / (grid.spacing * step));
        for s in 1..=reach.min(steps) {
            let by = [sign * d[0] * s, sign * d[1] * s, sign * d[2] * s];
            let Some(at) = grid.neighbour(cell, by) else {
                return false;
            };
            if solid[grid.in_cell(at)] {
                return true;
            }
        }
        false
    };
    let mut out = vec![false; grid.cells()];
    out.par_chunks_mut(grid.slab())
        .enumerate()
        .for_each(|(i, slab)| {
            if cancel.is_cancelled() {
                return;
            }
            for j in 0..grid.dims[1] {
                for k in 0..grid.dims[2] {
                    let at = grid.index(i, j, k);
                    if solid[at] || !reachable[at] {
                        continue;
                    }
                    let cell = [i, j, k];
                    let n = DIRECTIONS
                        .iter()
                        .filter(|d| hits(cell, **d, 1) && hits(cell, **d, -1))
                        .count();
                    if n >= options.min_directions {
                        slab[j * grid.dims[2] + k] = true;
                    }
                }
            }
        });
    cancel.check()?;
    Ok(out)
}

/// The 26-connected groups of marked cells with at least `min_points` cells, in scan order.
pub(super) fn components(
    grid: &Grid,
    marked: &[bool],
    min_points: usize,
    cancel: &CancellationToken,
) -> Result<Vec<Vec<[usize; 3]>>, Cancelled> {
    let mut seen = vec![false; grid.cells()];
    let mut groups: Vec<Vec<[usize; 3]>> = Vec::new();
    for i in 0..grid.dims[0] {
        cancel.check()?;
        for j in 0..grid.dims[1] {
            for k in 0..grid.dims[2] {
                let at = grid.index(i, j, k);
                if !marked[at] || seen[at] {
                    continue;
                }
                let mut stack = vec![[i, j, k]];
                let mut members = Vec::new();
                seen[at] = true;
                while let Some(cell) = stack.pop() {
                    members.push(cell);
                    for da in -1_i32..=1 {
                        for db in -1_i32..=1 {
                            for dc in -1_i32..=1 {
                                let Some(next) = grid.neighbour(cell, [da, db, dc]) else {
                                    continue;
                                };
                                let n = grid.in_cell(next);
                                if marked[n] && !seen[n] {
                                    seen[n] = true;
                                    stack.push(next);
                                }
                            }
                        }
                    }
                }
                if members.len() >= min_points {
                    groups.push(members);
                }
            }
        }
    }
    Ok(groups)
}
