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
use occlude_simd::NeighbourTile;

#[path = "stream/occlude_simd.rs"]
mod occlude_simd;

const BLOCK_ATOMS: usize = 64;
/// Smallest slab the memory fallback will shrink to before it refuses.
const MIN_SLAB_ATOMS: usize = 256;
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
        crate::neighbourhood::validate_positions(positions)?;
        let cutoff = self.cutoff;
        match CellList::build_in(
            positions,
            &self.targets,
            cutoff,
            periodic.copied(),
            CellGridOptions::default(),
            &self.context,
        ) {
            Ok(index) => self.sample_blocks(
                positions,
                periodic,
                &index,
                positions.len(),
                &|ordinal| ordinal,
                &mut emit,
            ),
            // The whole-structure index does not fit the budget. Degrade to
            // spatial slabs rather than refuse: peak memory then follows the slab
            // size, not the structure size, and areas stay bit-identical because
            // each atom still sees every neighbour within reach.
            Err(SpatialError::Memory(error)) if periodic.is_none() => {
                self.sample_slabs(positions, error, &mut emit)
            }
            Err(error) => Err(error.into()),
        }
    }

    /// Samples `count` atoms, named by `atom_at`, against one index and emits
    /// `(ordinal, area)` in ordinal order.
    fn sample_blocks(
        &self,
        positions: &[[f32; 3]],
        periodic: Option<&PeriodicBox>,
        index: &CellList<'_>,
        count: usize,
        atom_at: &(dyn Fn(usize) -> usize + Sync),
        emit: &mut impl FnMut(usize, f64) -> Result<(), SasaError>,
    ) -> Result<(), SasaError> {
        let context = &self.context;
        let sample_count = self.directions.len();
        let bytes_per_block = sample_count * size_of::<[f64; 3]>()
            + sample_count.div_ceil(64) * size_of::<u64>()
            + BLOCK_ATOMS * size_of::<f64>()
            + size_of::<[Neighbor; NEIGHBOR_TILE]>();
        let geometry = Sampling {
            positions,
            radii: self.radii,
            probe: self.probe,
            directions: &self.directions,
            cutoff: self.cutoff,
            index,
            periodic,
        };
        let mut next = 0;
        try_for_each_block_in(
            BlockPlan::new(count, BLOCK_ATOMS),
            context,
            bytes_per_block,
            |_, range| {
                let mut points = vec![[0.0; 3]; sample_count];
                let mut covered = vec![0_u64; sample_count.div_ceil(64)];
                let mut areas = Vec::with_capacity(range.len());
                for ordinal in range {
                    if context.cancellation().is_cancelled() {
                        return Err(SasaError::from(SpatialError::Cancelled));
                    }
                    areas.push(geometry.area(atom_at(ordinal), &mut points, &mut covered)?);
                }
                Ok(areas)
            },
            |areas| {
                for area in areas {
                    emit(next, area)?;
                    next += 1;
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

    /// Open-boundary fallback: sweep atoms in coordinate order along the longest
    /// axis, one slab of core atoms at a time. A slab's index holds the core plus
    /// a halo of one cutoff on each side, and the slab halves whenever its index
    /// is still refused, down to [`MIN_SLAB_ATOMS`]; `refusal` is returned when
    /// even that does not fit. Retained: one order array and one area per atom.
    fn sample_slabs(
        &self,
        positions: &[[f32; 3]],
        refusal: molframe_core::MemoryBudgetError,
        emit: &mut impl FnMut(usize, f64) -> Result<(), SasaError>,
    ) -> Result<(), SasaError> {
        let context = &self.context;
        let atoms = positions.len();
        let retained = atoms.saturating_mul(size_of::<u32>() + size_of::<f64>());
        let _reservation = context
            .try_reserve(retained)
            .map_err(|_| SasaError::from(SpatialError::Memory(refusal)))?;
        let axis = longest_axis(positions);
        let mut order: Vec<u32> = (0..atoms)
            .map(|atom| u32::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded))
            .collect::<Result<_, _>>()?;
        order.sort_unstable_by(|&a, &b| {
            positions[a as usize][axis]
                .total_cmp(&positions[b as usize][axis])
                .then(a.cmp(&b))
        });
        let coordinate = |atom: u32| f64::from(positions[atom as usize][axis]);
        let reach = f64::from(self.cutoff);
        let mut areas = vec![0.0_f64; atoms];
        let mut core = atoms.div_ceil(2).max(1);
        let mut start = 0;
        while start < atoms {
            if context.cancellation().is_cancelled() {
                return Err(SpatialError::Cancelled.into());
            }
            let end = (start + core).min(atoms);
            let low = coordinate(order[start]) - reach;
            let high = coordinate(order[end - 1]) + reach;
            let first = order.partition_point(|&atom| coordinate(atom) < low);
            let last = order.partition_point(|&atom| coordinate(atom) <= high);
            let index = match CellList::build_in(
                positions,
                &order[first..last],
                self.cutoff,
                None,
                CellGridOptions::default(),
                context,
            ) {
                Ok(index) => index,
                Err(SpatialError::Memory(_)) if core > MIN_SLAB_ATOMS => {
                    core = (core / 2).max(MIN_SLAB_ATOMS);
                    continue;
                }
                Err(SpatialError::Memory(_)) => {
                    return Err(SpatialError::Memory(refusal).into());
                }
                Err(error) => return Err(error.into()),
            };
            let slab = &order[start..end];
            self.sample_blocks(
                positions,
                None,
                &index,
                slab.len(),
                &|ordinal| slab[ordinal] as usize,
                &mut |ordinal, area| {
                    areas[slab[ordinal] as usize] = area;
                    Ok(())
                },
            )?;
            start = end;
        }
        for (atom, area) in areas.into_iter().enumerate() {
            emit(atom, area)?;
        }
        Ok(())
    }
}

/// The coordinate axis with the largest extent; ties go to the lowest axis.
fn longest_axis(positions: &[[f32; 3]]) -> usize {
    let mut low = [f32::INFINITY; 3];
    let mut high = [f32::NEG_INFINITY; 3];
    for position in positions {
        for axis in 0..3 {
            low[axis] = low[axis].min(position[axis]);
            high[axis] = high[axis].max(position[axis]);
        }
    }
    let extent = |axis: usize| high[axis] - low[axis];
    let mut best = 0;
    for axis in 1..3 {
        if extent(axis) > extent(best) {
            best = axis;
        }
    }
    best
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
        // Periodic sampling measures minimum-image distances and keeps the
        // scalar test; the open-boundary case uses the lane-parallel tile.
        let mut neighbors = [([0.0; 3], 0.0); NEIGHBOR_TILE];
        let mut neighbors_len = 0;
        let mut tile = NeighbourTile::<NEIGHBOR_TILE>::new();
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
                if self.periodic.is_some() {
                    neighbors[neighbors_len] = (other, other_radius * other_radius);
                    neighbors_len += 1;
                    if neighbors_len == NEIGHBOR_TILE {
                        hidden += self.occlude(points, covered, &neighbors);
                        neighbors_len = 0;
                    }
                } else {
                    tile.push(other, other_radius * other_radius);
                    if tile.is_full() {
                        hidden += tile.occlude(points, covered);
                        tile.clear();
                    }
                }
            })?;
        let sphere = 4.0 * core::f64::consts::PI;
        if !initialized {
            return Ok(sphere * radius * radius);
        }
        if neighbors_len != 0 {
            hidden += self.occlude(points, covered, &neighbors[..neighbors_len]);
        }
        if tile.len() != 0 {
            hidden += tile.occlude(points, covered);
        }
        Ok(f64::from(samples - hidden) * (sphere / f64::from(samples)) * radius * radius)
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
