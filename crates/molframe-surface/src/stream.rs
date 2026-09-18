//! Bounded Shrake–Rupley sampling without pair-proportional adjacency.
//!
//! A shared cell index supplies candidate neighbours. Each admitted block owns
//! one reusable sphere and occlusion bitmap, plus a bounded result tile. Runtime
//! is O(atoms × samples × local density); retained storage is O(atoms + workers
//! × samples), even when every atom occupies one dense cell.

use crate::accessible_area::{SasaError, point_on_sphere};
use crate::neighbourhood::{point_inside_sphere, squared_distance};
use crate::sampling_plan::SasaSampler;
use molframe_core::{
    ExecutionContext,
    parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in},
};
use molframe_spatial::{CellGridOptions, CellList, PeriodicBox, SpatialError};

const BLOCK_ATOMS: usize = 64;
const NEIGHBOR_TILE: usize = 64;
type Neighbor = ([f64; 3], f64);

/// Emits sampled atom areas in canonical atom order without storing adjacency.
///
/// Every allocation is reserved before construction; callbacks receive scalar
/// results. Sink-owned output needs its own reservation. Numerical tests
/// and sampling directions match the deterministic resident CPU reference.
/// Periodic sampling uses minimum-image distances to other atoms; self images
/// are excluded, matching the spatial pair convention.
///
/// # Errors
///
/// Returns invalid input, memory, cancellation, worker or consumer errors.
pub fn visit_shrake_rupley(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    samples: u16,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
    emit: impl FnMut(usize, f64) -> Result<(), SasaError>,
) -> Result<(), SasaError> {
    SasaSampler::new(radii, probe, samples, context)?.visit(positions, periodic, emit)
}

impl SasaSampler<'_> {
    /// Emits one frame's areas using cached radii validation, indices and directions.
    ///
    /// # Errors
    ///
    /// Returns dimension, geometry, resource, cancellation, worker or sink errors.
    pub fn visit(
        &self,
        positions: &[[f32; 3]],
        periodic: Option<&PeriodicBox>,
        mut emit: impl FnMut(usize, f64) -> Result<(), SasaError>,
    ) -> Result<(), SasaError> {
        if positions.len() != self.radii.len() {
            return Err(SasaError::LengthMismatch {
                positions: positions.len(),
                radii: self.radii.len(),
            });
        }
        let context = &self.context;
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled.into());
        }
        if positions.is_empty() {
            return Ok(());
        }
        let sample_count = self.directions.len();
        let cutoff = self.cutoff;
        let index = CellList::build_in(
            positions,
            &self.targets,
            cutoff,
            periodic.copied(),
            CellGridOptions::default(),
            context,
        )?;
        let bytes_per_block = sample_count * size_of::<[f64; 3]>()
            + sample_count.div_ceil(64) * size_of::<u64>()
            + BLOCK_ATOMS * size_of::<f64>()
            + size_of::<[Neighbor; NEIGHBOR_TILE]>();
        let geometry = Sampling {
            positions,
            radii: self.radii,
            probe: self.probe,
            directions: &self.directions,
            cutoff,
            index: &index,
            periodic,
        };
        let mut next_atom = 0;
        try_for_each_block_in(
            BlockPlan::new(positions.len(), BLOCK_ATOMS),
            context,
            bytes_per_block,
            |_, range| {
                let mut points = vec![[0.0; 3]; sample_count];
                let mut covered = vec![0_u64; sample_count.div_ceil(64)];
                let mut areas = Vec::with_capacity(range.len());
                for atom in range {
                    if context.cancellation().is_cancelled() {
                        return Err(SasaError::from(SpatialError::Cancelled));
                    }
                    areas.push(geometry.area(atom, &mut points, &mut covered)?);
                }
                Ok(areas)
            },
            |areas| {
                for area in areas {
                    emit(next_atom, area)?;
                    next_atom += 1;
                }
                Ok(())
            },
        )
        .map_err(|error| match error {
            BlockExecutionError::Memory(error) => SasaError::from(SpatialError::Memory(error)),
            BlockExecutionError::Cancelled => SasaError::from(SpatialError::Cancelled),
            BlockExecutionError::Worker(_) => SasaError::WorkerPanicked,
            BlockExecutionError::Operation(error) => error,
        })
    }
}

struct Sampling<'a> {
    positions: &'a [[f32; 3]],
    radii: &'a [f32],
    probe: f64,
    directions: &'a [[f64; 3]],
    cutoff: f32,
    index: &'a CellList<'a>,
    periodic: Option<&'a PeriodicBox>,
}

impl Sampling<'_> {
    // Test each uncovered point against a bounded tile, stopping at its first
    // occluder. This avoids scanning the bitmap once per neighbor while keeping
    // storage independent of the degree of a dense atom.
    fn occlude(&self, points: &[[f64; 3]], covered: &mut [u64], neighbors: &[Neighbor]) -> u16 {
        let mut hidden = 0;
        for (sample, &point) in points.iter().enumerate() {
            let mask = 1_u64 << (sample % 64);
            let word = &mut covered[sample / 64];
            if *word & mask == 0
                && neighbors.iter().any(|&(centre, radius_squared)| {
                    self.point_inside(point, centre, radius_squared)
                })
            {
                *word |= mask;
                hidden += 1;
            }
        }
        hidden
    }

    fn distance_squared(&self, left: [f64; 3], right: [f64; 3]) -> f64 {
        match self.periodic {
            Some(periodic) => squared_distance([0.0; 3], periodic.displacement_f64(left, right)),
            None => squared_distance(left, right),
        }
    }

    fn point_inside(&self, point: [f64; 3], centre: [f64; 3], radius_squared: f64) -> bool {
        match self.periodic {
            Some(_) => self.distance_squared(point, centre) < radius_squared,
            None => point_inside_sphere(point, centre, radius_squared),
        }
    }

    fn area(
        &self,
        atom: usize,
        points: &mut [[f64; 3]],
        covered: &mut [u64],
    ) -> Result<f64, SasaError> {
        let radius = f64::from(self.radii[atom]) + self.probe;
        if radius <= 0.0 {
            return Ok(0.0);
        }
        let centre = self.positions[atom].map(f64::from);
        let mut initialized = false;
        let mut hidden = 0_u16;
        let mut neighbors = [([0.0; 3], 0.0); NEIGHBOR_TILE];
        let mut neighbors_len = 0;
        let samples =
            u16::try_from(points.len()).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let atom = u32::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        self.index
            .for_each_candidate(&[atom], self.cutoff, |_, target, _| {
                if hidden == samples {
                    return;
                }
                let target = target as usize;
                let other = self.positions[target].map(f64::from);
                let other_radius = f64::from(self.radii[target]) + self.probe;
                let reach = radius + other_radius;
                if self.distance_squared(centre, other) >= reach * reach {
                    return;
                }
                if !initialized {
                    for (point, &direction) in points.iter_mut().zip(self.directions) {
                        *point = point_on_sphere(centre, radius, direction);
                    }
                    covered.fill(0);
                    initialized = true;
                }
                neighbors[neighbors_len] = (other, other_radius * other_radius);
                neighbors_len += 1;
                if neighbors_len == NEIGHBOR_TILE {
                    hidden += self.occlude(points, covered, &neighbors);
                    neighbors_len = 0;
                }
            })?;
        let sphere = 4.0 * core::f64::consts::PI;
        if !initialized {
            return Ok(sphere * radius * radius);
        }
        if neighbors_len != 0 {
            hidden += self.occlude(points, covered, &neighbors[..neighbors_len]);
        }
        Ok(f64::from(samples - hidden) * (sphere / f64::from(samples)) * radius * radius)
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
