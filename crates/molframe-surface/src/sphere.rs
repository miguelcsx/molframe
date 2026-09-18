//! A deterministic near-uniform set of directions on the unit sphere.
//!
//! Sampling accessibility needs points that are spread evenly and, above all,
//! the *same* every run — a random tessellation would make one structure's area
//! depend on a seed. The Fibonacci (golden-spiral) lattice places the `k`-th of
//! `n` points at height `1 - (2k+1)/n` and winds it around by the golden angle,
//! which gives a spacing close to optimal with no iteration and no lookup table.
//!
//! Cost is `O(n)` time and one buffer of `n` directions.

use core::f64::consts::PI;

/// Returns `count` unit vectors spread near-uniformly over the sphere.
///
/// The points are centred within their latitude bands so that none lands
/// exactly on a pole, and the sequence is fixed, so two runs with the same
/// `count` return byte-identical directions.
///
/// Runs in `O(count)` time.
#[must_use]
pub fn fibonacci_sphere(count: u16) -> Vec<[f64; 3]> {
    let mut directions = Vec::with_capacity(count as usize);

    for_each_fibonacci(count, |direction| {
        directions.push(direction);
    });

    directions
}

/// Streams deterministic Fibonacci-sphere directions to `visit`.
///
/// This provides the same sequence as [`fibonacci_sphere`] without allocating a
/// direction buffer, allowing per-atom density sampling to use `O(1)` temporary
/// space.
pub(crate) fn for_each_fibonacci(count: u16, mut visit: impl FnMut([f64; 3])) {
    if count == 0 {
        return;
    }

    let total = f64::from(count);
    let golden_angle = PI * (3.0 - 5.0_f64.sqrt());

    for index in 0..count {
        visit(fibonacci_direction(index, total, golden_angle));
    }
}

/// Computes one deterministic Fibonacci-lattice direction.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn fibonacci_direction(index: u16, total: f64, golden_angle: f64) -> [f64; 3] {
    let position = f64::from(index);
    let height = 1.0 - (2.0 * position + 1.0) / total;
    let radius = (1.0 - height * height).max(0.0).sqrt();
    let angle = golden_angle * position;

    [radius * angle.cos(), height, radius * angle.sin()]
}

#[cfg(test)]
#[path = "sphere_tests.rs"]
mod tests;
