//! Lane-parallel occlusion of sample points by a tile of neighbour spheres.
//!
//! A sample point is hidden when it lies strictly inside any neighbour's
//! expanded sphere. The neighbours of one tile sit in structure-of-arrays form,
//! so one sample point is tested against several of them per instruction. The
//! squared distance is formed exactly as the scalar test forms it — three
//! differences, three products, two additions left to right, no fused
//! multiply-add — so the answer is bit-identical to the scalar loop.

#[cfg(target_arch = "x86_64")]
use wide::f64x4 as F64xN;
#[cfg(target_arch = "x86_64")]
const LANES: usize = 4;
#[cfg(not(target_arch = "x86_64"))]
use wide::f64x2 as F64xN;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 2;

/// Neighbour spheres held one coordinate per array.
///
/// Slots past `len` up to the next lane boundary hold a zero squared radius,
/// which no squared distance is strictly below, so a padded lane never hides
/// a point.
pub(super) struct NeighbourTile<const N: usize> {
    x: [f64; N],
    y: [f64; N],
    z: [f64; N],
    radius_squared: [f64; N],
    len: usize,
}

impl<const N: usize> NeighbourTile<N> {
    pub(super) const fn new() -> Self {
        Self {
            x: [0.0; N],
            y: [0.0; N],
            z: [0.0; N],
            radius_squared: [0.0; N],
            len: 0,
        }
    }

    pub(super) const fn len(&self) -> usize {
        self.len
    }

    pub(super) const fn is_full(&self) -> bool {
        self.len == N
    }

    /// Adds one sphere; the caller flushes the tile before it overflows.
    pub(super) fn push(&mut self, centre: [f64; 3], radius_squared: f64) {
        let slot = self.len;
        self.x[slot] = centre[0];
        self.y[slot] = centre[1];
        self.z[slot] = centre[2];
        self.radius_squared[slot] = radius_squared;
        self.len += 1;
    }

    /// Empties the tile, restoring the never-hiding padding.
    pub(super) fn clear(&mut self) {
        self.radius_squared[..self.len].fill(0.0);
        self.len = 0;
    }

    /// Marks every uncovered sample the tile hides and returns how many it hid.
    pub(super) fn occlude(&self, points: &[[f64; 3]], covered: &mut [u64]) -> u16 {
        let padded = self.len.div_ceil(LANES) * LANES;
        let mut hidden = 0;
        for (sample, point) in points.iter().enumerate() {
            let mask = 1_u64 << (sample % 64);
            let word = &mut covered[sample / 64];
            if *word & mask == 0 && self.hides(*point, padded) {
                *word |= mask;
                hidden += 1;
            }
        }
        hidden
    }

    fn hides(&self, point: [f64; 3], padded: usize) -> bool {
        let px = F64xN::splat(point[0]);
        let py = F64xN::splat(point[1]);
        let pz = F64xN::splat(point[2]);
        let mut start = 0;
        while start < padded {
            let lanes = |values: &[f64; N]| {
                F64xN::from(core::array::from_fn::<f64, LANES, _>(|lane| {
                    values[start + lane]
                }))
            };
            let dx = px - lanes(&self.x);
            let dy = py - lanes(&self.y);
            let dz = pz - lanes(&self.z);
            let distance_squared = dx * dx + dy * dy + dz * dz;
            if distance_squared.simd_lt(lanes(&self.radius_squared)).any() {
                return true;
            }
            start += LANES;
        }
        false
    }
}

#[cfg(test)]
#[path = "occlude_simd_tests.rs"]
mod tests;
