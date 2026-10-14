//! Distribution probe that keeps `auto` away from cell lists on clustered data.

use crate::brute::finite;
use molframe_core::numeric::{f32_to_i64, usize_to_f64};

/// Largest sample the probe inspects, so its cost is independent of system size.
const SAMPLE: usize = 65_536;
/// Mean points sharing a cell with a point, beyond which `auto` prefers a k-d tree.
///
/// A cell list tests each point against every point in 27 cells, so work grows
/// with the mean occupancy seen by points, while the k-d tree adapts to local
/// density and scans leaves in parallel. Measured: on uniform clouds the tree
/// matches the cell list in CPU time near 50 points per cell and is faster in
/// wall time everywhere; on a clustered cloud with a dense core (about 12
/// here) it is 1.5x cheaper in CPU and 4x in wall time. Uniform clouds below
/// the limit stay on the cell list.
const CROWDING_LIMIT: f64 = 8.0;

/// Whether cutoff-sized cells hold so many points on average that a k-d tree
/// does less work than a cell list.
///
/// Estimates the occupancy seen by an average point, `sum(c^2) / n`, from an
/// evenly strided sample. Runtime is `O(S log S)` for sample size `S`.
pub(super) fn is_crowded(positions: &[[f32; 3]], atoms: &[u32], cutoff: f32) -> bool {
    if atoms.len() < 2 * 1_024 || !cutoff.is_finite() || cutoff <= 0.0 {
        return false;
    }
    let stride = atoms.len().div_ceil(SAMPLE);
    let mut keys: Vec<[i64; 3]> = atoms
        .iter()
        .step_by(stride)
        .filter_map(|&atom| positions.get(atom as usize).copied())
        .filter(|&point| finite(point))
        .map(|point| point.map(|value| f32_to_i64((value / cutoff).floor())))
        .collect();
    if keys.is_empty() {
        return false;
    }
    keys.sort_unstable();
    let mut sum_squares = 0.0_f64;
    let mut run = 0.0_f64;
    let mut previous = keys[0];
    for key in &keys {
        if *key == previous {
            run += 1.0;
        } else {
            sum_squares += run * run;
            run = 1.0;
            previous = *key;
        }
    }
    sum_squares += run * run;
    let sampled = usize_to_f64(keys.len());
    let rate = 1.0 / usize_to_f64(stride);
    // Unbiased estimate of sum(c^2) / n: a sampled count c_s has
    // E[c_s^2 - c_s] = rate^2 * c * (c - 1), which removes the sampling noise
    // that would otherwise add `1 / rate` to every estimate.
    (sum_squares - sampled) / (sampled * rate) + 1.0 > CROWDING_LIMIT
}
