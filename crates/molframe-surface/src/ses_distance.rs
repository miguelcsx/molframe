//! Exact linear-time Euclidean distance transform for the SES voxel field.
//!
//! The three-dimensional transform is separated into one-dimensional lower
//! envelopes along X, Y and Z. Each cell is read and written once per axis;
//! scratch storage is reused and depends only on the longest grid dimension.

use crate::SasaError;
use crate::cavity::{EXTERIOR, Grid};
use crate::numeric::{f64_to_f32, usize_to_f64, usize_to_u32};
use crate::workspace::empty_with_capacity;

/// Computes exact Euclidean distance from every cell to the exterior set.
pub(super) fn distance_from_exterior(grid: &Grid, state: Vec<u8>) -> Result<Vec<f32>, SasaError> {
    let mut distance = initial_distances(&state)?;
    drop(state);

    let maximum_line = grid
        .dims
        .into_iter()
        .max()
        .ok_or(SasaError::GridDimensionsOverflow)?;
    let mut scratch = DistanceScratch::new(maximum_line)?;
    transform_x(grid, &mut distance, &mut scratch);
    transform_y(grid, &mut distance, &mut scratch);
    transform_z(grid, &mut distance, &mut scratch);
    drop(scratch);

    let step = f64_to_f32(grid.step);
    for value in &mut distance {
        *value = value.sqrt() * step;
    }
    Ok(distance)
}

fn initial_distances(state: &[u8]) -> Result<Vec<f32>, SasaError> {
    let mut distance = empty_with_capacity(state.len())?;
    for &cell in state {
        distance.push(if cell == EXTERIOR { 0.0 } else { f32::INFINITY });
    }
    Ok(distance)
}

/// Reusable lower-envelope storage for one grid line.
struct DistanceScratch {
    sites: Vec<u32>,
    costs: Vec<f32>,
    starts: Vec<f64>,
}

impl DistanceScratch {
    fn new(maximum_line: usize) -> Result<Self, SasaError> {
        Ok(Self {
            sites: empty_with_capacity(maximum_line)?,
            costs: empty_with_capacity(maximum_line)?,
            starts: empty_with_capacity(maximum_line)?,
        })
    }

    /// Applies the squared-distance transform to one strided line in place.
    fn transform_line(&mut self, field: &mut [f32], start: usize, stride: usize, length: usize) {
        self.sites.clear();
        self.costs.clear();
        self.starts.clear();

        for coordinate in 0..length {
            let value = field[start + coordinate * stride];
            if value.is_finite() {
                self.add_site(coordinate, value);
            }
        }
        if self.sites.is_empty() {
            return;
        }

        let mut envelope = 0usize;
        for coordinate in 0..length {
            let coordinate_f64 = usize_to_f64(coordinate);
            while envelope + 1 < self.sites.len() && self.starts[envelope + 1] < coordinate_f64 {
                envelope += 1;
            }
            let site = u32_to_usize(self.sites[envelope]);
            let delta = coordinate_f64 - usize_to_f64(site);
            field[start + coordinate * stride] =
                f64_to_f32(delta * delta + f64::from(self.costs[envelope]));
        }
    }

    /// Adds one parabola while maintaining the lower envelope.
    fn add_site(&mut self, coordinate: usize, cost: f32) {
        loop {
            let Some(&last_site) = self.sites.last() else {
                self.sites.push(usize_to_u32(coordinate));
                self.costs.push(cost);
                self.starts.push(f64::NEG_INFINITY);
                return;
            };
            let last = self.sites.len() - 1;
            let crossing =
                intersection(u32_to_usize(last_site), self.costs[last], coordinate, cost);
            if crossing > self.starts[last] {
                self.sites.push(usize_to_u32(coordinate));
                self.costs.push(cost);
                self.starts.push(crossing);
                return;
            }
            self.sites.pop();
            self.costs.pop();
            self.starts.pop();
        }
    }
}

fn u32_to_usize(value: u32) -> usize {
    match usize::try_from(value) {
        Ok(converted) => converted,
        Err(_) => usize::MAX,
    }
}

fn intersection(first: usize, first_cost: f32, second: usize, second_cost: f32) -> f64 {
    let first = usize_to_f64(first);
    let second = usize_to_f64(second);
    (f64::from(second_cost) + second * second - f64::from(first_cost) - first * first)
        / (2.0 * (second - first))
}

fn transform_x(grid: &Grid, field: &mut [f32], scratch: &mut DistanceScratch) {
    for z in 0..grid.dims[2] {
        for y in 0..grid.dims[1] {
            scratch.transform_line(field, grid.index(0, y, z), 1, grid.dims[0]);
        }
    }
}

fn transform_y(grid: &Grid, field: &mut [f32], scratch: &mut DistanceScratch) {
    let stride = grid.dims[0];
    for z in 0..grid.dims[2] {
        for x in 0..grid.dims[0] {
            scratch.transform_line(field, grid.index(x, 0, z), stride, grid.dims[1]);
        }
    }
}

fn transform_z(grid: &Grid, field: &mut [f32], scratch: &mut DistanceScratch) {
    let stride = grid.dims[0] * grid.dims[1];
    for y in 0..grid.dims[1] {
        for x in 0..grid.dims[0] {
            scratch.transform_line(field, grid.index(x, y, 0), stride, grid.dims[2]);
        }
    }
}

#[cfg(test)]
#[path = "ses_distance_tests.rs"]
mod tests;
