//! A provable upper bound on the TM-score from intra-structure distances.
//!
//! A rigid motion preserves distances inside each structure, so for any two
//! points the superposed residuals obey `d_i + d_j >= delta_ij`, where
//! `delta_ij` is the difference between the point pair's distance in the model
//! and in the reference. Points whose pairwise `delta` is at least `D` cannot
//! both lie closer than `D / 2` to their targets, which caps their summed TM
//! kernel. Splitting the points into such groups gives a bound on the best
//! score of any superposition, and therefore on [`crate::tm_score`], without
//! fitting anything.

use crate::numeric::usize_to_f64;
use crate::superposed::tm_scale;

/// Most points the bound inspects; the rest each count as a full kernel.
const MAX_LANDMARKS: usize = 512;

/// Group separations tried, as multiples of the distance scale.
const SEPARATIONS: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 6.0];

/// Absolute slack on the normalised bound for coordinate rounding.
const SLACK: f64 = 1e-4;

const NONE: u32 = u32::MAX;

/// Reusable storage for [`TmBound::bound`].
#[derive(Default)]
pub(crate) struct TmBound {
    delta: Vec<f64>,
    head: Vec<u32>,
    next: Vec<u32>,
}

impl TmBound {
    /// Bytes of storage one bound needs for structures of `length` points.
    pub(crate) fn bytes(length: usize) -> usize {
        let points = length.min(MAX_LANDMARKS);
        points * points * size_of::<f64>() + points * 2 * size_of::<u32>()
    }

    /// An upper bound on the TM-score of `model` against `reference`, which
    /// must have equal lengths of at least one point.
    ///
    /// No superposition of the pair scores above the returned value.
    pub(crate) fn bound(&mut self, model: &[[f32; 3]], reference: &[[f32; 3]]) -> f64 {
        let length = model.len();
        let points = length.min(MAX_LANDMARKS);
        let picked = |slot: usize| slot * length / points;
        self.delta.clear();
        self.delta.resize(points * points, 0.0);
        for row in 0..points {
            for column in row + 1..points {
                let (i, j) = (picked(row), picked(column));
                let value = (gap(model[i], model[j]) - gap(reference[i], reference[j])).abs();
                self.delta[row * points + column] = value;
                self.delta[column * points + row] = value;
            }
        }
        let scale = tm_scale(length);
        let mut kernel_sum = usize_to_f64(points);
        for separation in SEPARATIONS {
            let distance = separation * scale;
            let groups = usize_to_f64(self.group(points, distance));
            let half = distance / 2.0 / scale;
            let kernel = 1.0 / (1.0 + half * half);
            kernel_sum = kernel_sum.min(groups + (usize_to_f64(points) - groups) * kernel);
        }
        let outside = usize_to_f64(length - points);
        ((outside + kernel_sum) / usize_to_f64(length) + SLACK).min(1.0)
    }

    /// Greedily covers the points with groups whose members all differ by at
    /// least `distance`, returning the group count.
    fn group(&mut self, points: usize, distance: f64) -> usize {
        self.head.clear();
        self.head.resize(points, NONE);
        self.next.clear();
        self.next.resize(points, NONE);
        let mut groups = 0;
        for point in 0..points {
            let row = &self.delta[point * points..(point + 1) * points];
            let mut placed = false;
            for group in 0..groups {
                if members(&self.next, self.head[group]).all(|member| row[member] >= distance) {
                    self.next[point] = self.head[group];
                    self.head[group] = to_link(point);
                    placed = true;
                    break;
                }
            }
            if !placed {
                self.head[groups] = to_link(point);
                groups += 1;
            }
        }
        groups
    }
}

fn to_link(point: usize) -> u32 {
    u32::try_from(point).unwrap_or(NONE)
}

/// The members of one group, newest first.
fn members(next: &[u32], head: u32) -> impl Iterator<Item = usize> + '_ {
    let mut cursor = head;
    std::iter::from_fn(move || {
        let current = usize::try_from(cursor).ok().filter(|_| cursor != NONE)?;
        cursor = next[current];
        Some(current)
    })
}

fn gap(a: [f32; 3], b: [f32; 3]) -> f64 {
    molframe_geom::distance(a, b)
}

#[cfg(test)]
#[path = "tm_bound_tests.rs"]
pub(crate) mod tests;
