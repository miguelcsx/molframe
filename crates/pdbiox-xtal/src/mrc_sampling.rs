//! Coordinate-space interpolation over canonical density grids.

use super::{DensityMap, MapBoundary, linear, resolve_index};
use crate::CellTransform;
use crate::numeric::{f64_to_f32, f64_to_i64, usize_to_f64};

impl DensityMap {
    /// Returns the X-fastest scalar at one local grid index.
    #[must_use]
    pub fn value(&self, index: [usize; 3]) -> Option<f32> {
        if index
            .iter()
            .zip(self.dimensions)
            .any(|(index, dimension)| *index >= dimension)
        {
            return None;
        }
        self.values.get(linear(index, self.dimensions)).copied()
    }

    /// Trilinearly samples a Cartesian coordinate in ångströms.
    #[must_use]
    pub fn sample_cartesian(&self, position: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        self.sample_grid(self.cartesian_to_grid(position)?, boundary)
    }

    /// Tricubically samples a Cartesian coordinate with a Catmull-Rom kernel.
    #[must_use]
    pub fn sample_cartesian_cubic(&self, position: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        self.sample_grid_cubic(self.cartesian_to_grid(position)?, boundary)
    }

    /// Prepares the cell transform once so many samples share it.
    ///
    /// Building the transform costs several transcendental functions and a
    /// matrix inverse — an order of magnitude more than the interpolation it
    /// enables — so sampling a coordinate list one call at a time spends nearly
    /// all its time rebuilding a value that never changes.
    #[must_use]
    pub fn sampler(&self) -> Option<DensitySampler<'_>> {
        Some(DensitySampler {
            map: self,
            transform: CellTransform::new(&self.cell).ok()?,
        })
    }

    fn cartesian_to_grid(&self, position: [f64; 3]) -> Option<[f64; 3]> {
        let transform = CellTransform::new(&self.cell).ok()?;
        let has_origin = self.origin.iter().any(|value| value.abs() > f64::EPSILON);
        let fractional = if has_origin {
            transform.to_fractional([
                position[0] - self.origin[0],
                position[1] - self.origin[1],
                position[2] - self.origin[2],
            ])
        } else {
            transform.to_fractional(position)
        };
        let mut local = [0.0; 3];
        for axis in 0..3 {
            local[axis] = fractional[axis] * usize_to_f64(self.sampling[axis]);
            if !has_origin {
                local[axis] -= f64::from(self.starts[axis]);
            }
        }
        Some(local)
    }

    /// Trilinearly samples a local grid coordinate.
    #[must_use]
    pub fn sample_grid(&self, coordinate: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        if coordinate.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let base = coordinate.map(f64::floor);
        let fraction = [
            coordinate[0] - base[0],
            coordinate[1] - base[1],
            coordinate[2] - base[2],
        ];
        let mut result = 0.0_f64;
        for dz in 0..=1 {
            for dy in 0..=1 {
                for dx in 0..=1 {
                    let weight = axis_weight(fraction[0], dx)
                        * axis_weight(fraction[1], dy)
                        * axis_weight(fraction[2], dz);
                    if weight == 0.0 {
                        continue;
                    }
                    let index = resolve_index(
                        [
                            f64_to_i64(base[0]) + dx,
                            f64_to_i64(base[1]) + dy,
                            f64_to_i64(base[2]) + dz,
                        ],
                        self.dimensions,
                        boundary,
                    )?;
                    result += f64::from(self.values[linear(index, self.dimensions)]) * weight;
                }
            }
        }
        Some(f64_to_f32(result))
    }

    /// Tricubically samples a local grid coordinate with Catmull-Rom weights.
    #[must_use]
    pub fn sample_grid_cubic(&self, coordinate: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        if coordinate.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let base = coordinate.map(f64::floor);
        let weights = [
            cubic_weights(coordinate[0] - base[0]),
            cubic_weights(coordinate[1] - base[1]),
            cubic_weights(coordinate[2] - base[2]),
        ];
        let offsets = [-1_i64, 0, 1, 2];
        let mut result = 0.0_f64;
        for (z_slot, z_offset) in offsets.into_iter().enumerate() {
            for (y_slot, y_offset) in offsets.into_iter().enumerate() {
                for (x_slot, x_offset) in offsets.into_iter().enumerate() {
                    let weight = weights[0][x_slot] * weights[1][y_slot] * weights[2][z_slot];
                    if weight == 0.0 {
                        continue;
                    }
                    let index = resolve_index(
                        [
                            f64_to_i64(base[0]) + x_offset,
                            f64_to_i64(base[1]) + y_offset,
                            f64_to_i64(base[2]) + z_offset,
                        ],
                        self.dimensions,
                        boundary,
                    )?;
                    result += f64::from(self.values[linear(index, self.dimensions)]) * weight;
                }
            }
        }
        Some(f64_to_f32(result))
    }
}

fn axis_weight(fraction: f64, offset: i64) -> f64 {
    if offset == 0 {
        1.0 - fraction
    } else {
        fraction
    }
}

fn cubic_weights(fraction: f64) -> [f64; 4] {
    let squared = fraction * fraction;
    let cubed = squared * fraction;
    [
        -0.5 * fraction + squared - 0.5 * cubed,
        1.0 - 2.5 * squared + 1.5 * cubed,
        0.5 * fraction + 2.0 * squared - 1.5 * cubed,
        -0.5 * squared + 0.5 * cubed,
    ]
}

/// A density map with its cell transform already built.
///
/// Sampling through this pays for the transform once instead of once per
/// coordinate, and holds nothing proportional to the number of samples.
#[derive(Debug)]
pub struct DensitySampler<'a> {
    /// The map being sampled.
    map: &'a DensityMap,
    /// The prepared Cartesian-to-fractional conversion.
    transform: CellTransform,
}

impl DensitySampler<'_> {
    /// Trilinearly samples a Cartesian coordinate in ångströms.
    #[must_use]
    pub fn sample_cartesian(&self, position: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        self.map.sample_grid(self.to_grid(position), boundary)
    }

    /// Tricubically samples a Cartesian coordinate with a Catmull-Rom kernel.
    #[must_use]
    pub fn sample_cartesian_cubic(&self, position: [f64; 3], boundary: MapBoundary) -> Option<f32> {
        self.map.sample_grid_cubic(self.to_grid(position), boundary)
    }

    /// Converts a Cartesian coordinate to a local grid coordinate.
    ///
    /// Identical arithmetic to `DensityMap::cartesian_to_grid`, reading the
    /// prepared transform instead of constructing one.
    fn to_grid(&self, position: [f64; 3]) -> [f64; 3] {
        let map = self.map;
        let has_origin = map.origin.iter().any(|value| value.abs() > f64::EPSILON);
        let fractional = if has_origin {
            self.transform.to_fractional([
                position[0] - map.origin[0],
                position[1] - map.origin[1],
                position[2] - map.origin[2],
            ])
        } else {
            self.transform.to_fractional(position)
        };

        let mut local = [0.0; 3];
        for axis in 0..3 {
            local[axis] = fractional[axis] * usize_to_f64(map.sampling[axis]);
            if !has_origin {
                local[axis] -= f64::from(map.starts[axis]);
            }
        }
        local
    }
}
