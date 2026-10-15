//! A point set bucketed by cell, for asking how close the nearest of its points is.

use crate::numeric::f32_to_i64;
use std::collections::HashMap;

/// Cells a little wider than the radius, so rounding at a cell edge cannot hide a point.
const SLACK: f32 = 1.001;

/// Points bucketed by cell, with the box they fill, for proximity queries within one radius.
#[derive(Clone, Debug)]
pub struct PointProximity {
    lo: [f32; 3],
    hi: [f32; 3],
    edge: f32,
    cells: HashMap<[i64; 3], Vec<[f32; 3]>>,
}

impl PointProximity {
    /// Buckets `points` for queries about distances up to `radius`, ångström.
    ///
    /// Runs in `O(points)` time.
    #[must_use]
    pub fn new(points: &[[f32; 3]], radius: f32) -> Self {
        let edge = (radius * SLACK).max(f32::MIN_POSITIVE);
        let (mut lo, mut hi) = ([f32::MAX; 3], [f32::MIN; 3]);
        let mut cells: HashMap<[i64; 3], Vec<[f32; 3]>> = HashMap::new();
        for p in points {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
            cells.entry(Self::cell(*p, edge)).or_default().push(*p);
        }
        Self {
            lo,
            hi,
            edge,
            cells,
        }
    }

    fn cell(p: [f32; 3], edge: f32) -> [i64; 3] {
        p.map(|c| f32_to_i64((c / edge).floor()))
    }

    /// Whether `p` can be within `radius` of a point: inside the points' box, padded.
    #[must_use]
    pub fn may_reach(&self, p: [f32; 3], radius: f32) -> bool {
        let pad = radius * SLACK;
        (0..3).all(|k| p[k] >= self.lo[k] - pad && p[k] <= self.hi[k] + pad)
    }

    /// The distance from `p` to the nearest point in the cells around it, which holds every
    /// point within the radius given to [`new`](Self::new); `f32::MAX` when there is none.
    #[must_use]
    pub fn nearest(&self, p: [f32; 3]) -> f32 {
        let c = Self::cell(p, self.edge);
        let mut nearest = f32::MAX;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = self.cells.get(&[
                        c[0].saturating_add(dx),
                        c[1].saturating_add(dy),
                        c[2].saturating_add(dz),
                    ]) else {
                        continue;
                    };
                    for q in bucket {
                        let d =
                            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2))
                                .sqrt();
                        nearest = nearest.min(d);
                    }
                }
            }
        }
        nearest
    }
}
