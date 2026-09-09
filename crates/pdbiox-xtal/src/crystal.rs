//! Symmetry-aware crystal neighbour enumeration with bounded output.

use crate::SymmetrySet;
use crate::crystal_batch::visit_prepared_image_batches;
use crate::crystal_images::{CrystalImage, PreparedImages, SourcePosition, prepare_images};
use crate::numeric::f64_to_f32;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::{
    AtomIndex, Batch, BatchDemand, BatchLease, Code, Diagnostic, ExecutionContext,
    MemoryReservation, ModelIndex, Structure,
};
use pdbiox_spatial::{
    PairQuery, SpatialBackend, SpatialSearchOptions, for_each_pairs_within_unsorted,
};

/// Default ceiling on candidate atom images examined by a crystal search.
pub const DEFAULT_CRYSTAL_IMAGE_LIMIT: usize = 10_000_000;

/// One unique contact from an asymmetric-unit atom to a crystal image.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystalNeighbor {
    /// Atom anchored in the deposited asymmetric unit.
    pub source_atom: AtomIndex,
    /// Source atom transformed into the neighbouring image.
    pub image_atom: AtomIndex,
    /// Position of the explicit symmetry representative.
    pub operation: usize,
    /// Primitive lattice translation after applying the representative.
    pub lattice: [i32; 3],
    /// Squared Cartesian distance in ångström².
    pub distance_squared: f64,
}

/// One caller-bounded run of crystal neighbours.
#[derive(Debug)]
pub struct CrystalNeighborBatch {
    neighbors: Vec<CrystalNeighbor>,
}

/// Planning, work and batch controls for crystal-neighbour enumeration.
#[derive(Clone, Copy, Debug)]
pub struct CrystalNeighborOptions {
    /// Spatial search backend.
    pub backend: SpatialBackend,
    /// Maximum candidate images examined.
    pub candidate_limit: usize,
    /// Maximum rows and retained bytes per neighbour batch.
    pub demand: BatchDemand,
}

impl Default for CrystalNeighborOptions {
    fn default() -> Self {
        Self {
            backend: SpatialBackend::Auto,
            candidate_limit: DEFAULT_CRYSTAL_IMAGE_LIMIT,
            demand: BatchDemand::new(4096, 8_000_000),
        }
    }
}

impl CrystalNeighborBatch {
    /// Neighbours in deterministic generation order.
    #[must_use]
    pub fn neighbors(&self) -> &[CrystalNeighbor] {
        &self.neighbors
    }
}

impl Batch for CrystalNeighborBatch {
    fn rows(&self) -> usize {
        self.neighbors.len()
    }

    fn retained_bytes(&self) -> usize {
        self.neighbors
            .capacity()
            .saturating_mul(std::mem::size_of::<CrystalNeighbor>())
    }
}

/// Visits each unique neighbour without retaining the complete result.
///
/// Reverse-equivalent contacts are rejected during generation, so the path is
/// O(candidate pairs + output) and needs neither a global set nor a final sort.
///
/// # Errors
///
/// Returns diagnostics for invalid crystal data, resources, cancellation, or
/// spatial-search failure.
pub fn visit_crystal_neighbors(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    options: CrystalNeighborOptions,
    context: &ExecutionContext,
    mut visitor: impl FnMut(CrystalNeighbor) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let prepared = prepare_images(structure, symmetry, model, cutoff, context)?;
    let available = context
        .memory_budget()
        .bytes()
        .saturating_sub(context.reserved_bytes());
    let image_bytes = available
        .saturating_div(4)
        .max(std::mem::size_of::<CrystalImage>());
    visit_prepared_image_batches(
        &prepared,
        symmetry,
        options.candidate_limit,
        BatchDemand::new(4096, image_bytes),
        context,
        |images| {
            visit_image_neighbors(
                &prepared,
                images.batch().images(),
                symmetry,
                cutoff,
                options.backend,
                context,
                &mut visitor,
            )
        },
    )
}

/// Visits budget-charged neighbour batches.
///
/// # Errors
///
/// Returns a resource diagnostic when demand cannot retain one neighbour or
/// the context cannot reserve the next batch.
pub fn crystal_neighbor_batches(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    options: CrystalNeighborOptions,
    context: &ExecutionContext,
    mut visitor: impl FnMut(BatchLease<CrystalNeighborBatch>) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let mut accumulator = NeighborAccumulator::new(options.demand, context)?;
    visit_crystal_neighbors(
        structure,
        symmetry,
        model,
        cutoff,
        options,
        context,
        |neighbor| accumulator.push(neighbor, context, &mut visitor),
    )?;
    accumulator.finish(&mut visitor)
}

/// Explicitly materialises all neighbours under the shared memory budget.
///
/// # Errors
///
/// Returns before the output grows beyond the remaining context capacity.
pub fn collect_crystal_neighbors(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    options: CrystalNeighborOptions,
    context: &ExecutionContext,
) -> Result<Vec<CrystalNeighbor>, Diagnostic> {
    let mut output = Vec::new();
    let width = std::mem::size_of::<CrystalNeighbor>();
    let mut reservation = reserve(context, 0)?;
    visit_crystal_neighbors(
        structure,
        symmetry,
        model,
        cutoff,
        options,
        context,
        |neighbor| {
            if output.len() == output.capacity() {
                let additional = output.capacity().max(256);
                let additional_bytes = additional.checked_mul(width).ok_or_else(search_limit)?;
                reservation
                    .try_grow(additional_bytes)
                    .map_err(|error| resource_error(&error.to_string()))?;
                output
                    .try_reserve_exact(additional)
                    .map_err(|error| resource_error(&error.to_string()))?;
            }
            output.push(neighbor);
            Ok(())
        },
    )?;
    let _reservation = reservation;
    Ok(output)
}

fn visit_image_neighbors(
    prepared: &PreparedImages,
    images: &[CrystalImage],
    symmetry: &SymmetrySet,
    cutoff: f64,
    backend: SpatialBackend,
    context: &ExecutionContext,
    visitor: &mut impl FnMut(CrystalNeighbor) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let source_count = u32::try_from(prepared.sources.len()).map_err(|_| search_limit())?;
    let total = prepared
        .sources
        .len()
        .checked_add(images.len())
        .ok_or_else(search_limit)?;
    let total_count = u32::try_from(total).map_err(|_| search_limit())?;
    let bytes = total
        .checked_mul(std::mem::size_of::<[f32; 3]>())
        .ok_or_else(search_limit)?;
    let _reservation = reserve(context, bytes)?;
    let mut positions = Vec::with_capacity(total);
    positions.extend(prepared.sources.iter().map(|source| source.cartesian));
    positions.extend(images.iter().map(|image| image.cartesian));
    let left = AtomSelection::All(source_count);
    let right = AtomSelection::range(source_count..total_count);
    let mut failure = None;
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions: &positions,
            left: &left,
            right: &right,
            cutoff: f64_to_f32(cutoff),
            options: SpatialSearchOptions::with_backend(backend),
            periodic: None,
            context,
        },
        |pair| {
            if failure.is_some() {
                return;
            }
            match crystal_neighbor(pair, &prepared.sources, images, symmetry, source_count) {
                Ok(Some(neighbor)) => {
                    if let Err(error) = visitor(neighbor) {
                        failure = Some(error);
                    }
                }
                Ok(None) => {}
                Err(error) => failure = Some(error),
            }
        },
    )
    .map_err(pdbiox_spatial::SpatialError::into_diagnostic)?;
    match failure {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PairKey {
    source: u32,
    image: u32,
    operation: u32,
    lattice: [i32; 3],
}

impl PairKey {
    fn new(
        source: AtomIndex,
        image: AtomIndex,
        operation: usize,
        lattice: [i32; 3],
    ) -> Result<Self, Diagnostic> {
        Ok(Self {
            source: source.get(),
            image: image.get(),
            operation: u32::try_from(operation).map_err(|_| search_limit())?,
            lattice,
        })
    }
}

fn crystal_neighbor(
    pair: pdbiox_spatial::NeighborPair,
    sources: &[SourcePosition],
    images: &[CrystalImage],
    symmetry: &SymmetrySet,
    source_count: u32,
) -> Result<Option<CrystalNeighbor>, Diagnostic> {
    let source = sources.get(pair.first as usize).ok_or_else(invariant)?;
    let image_index = pair
        .second
        .checked_sub(source_count)
        .ok_or_else(invariant)?;
    let image = images.get(image_index as usize).ok_or_else(invariant)?;
    let distance_squared = f64::from(pair.distance_squared);
    if source.atom == image.atom && distance_squared <= f64::EPSILON {
        return Ok(None);
    }
    let forward = PairKey::new(source.atom, image.atom, image.operation, image.lattice)?;
    let (inverse_operation, inverse_lattice) =
        symmetry.inverse_image(image.operation, image.lattice)?;
    let reverse = PairKey::new(image.atom, source.atom, inverse_operation, inverse_lattice)?;
    if forward > reverse {
        return Ok(None);
    }
    Ok(Some(CrystalNeighbor {
        source_atom: source.atom,
        image_atom: image.atom,
        operation: image.operation,
        lattice: image.lattice,
        distance_squared,
    }))
}

struct NeighborAccumulator {
    values: Vec<CrystalNeighbor>,
    reservation: Option<MemoryReservation>,
    capacity: usize,
}

impl NeighborAccumulator {
    fn new(demand: BatchDemand, context: &ExecutionContext) -> Result<Self, Diagnostic> {
        let width = std::mem::size_of::<CrystalNeighbor>();
        let capacity = demand.max_rows.min(demand.max_bytes / width);
        if capacity == 0 {
            return Err(resource_error(
                "batch demand cannot retain one crystal neighbor",
            ));
        }
        let bytes = capacity.checked_mul(width).ok_or_else(search_limit)?;
        Ok(Self {
            values: Vec::with_capacity(capacity),
            reservation: Some(reserve(context, bytes)?),
            capacity,
        })
    }

    fn push(
        &mut self,
        value: CrystalNeighbor,
        context: &ExecutionContext,
        visitor: &mut impl FnMut(BatchLease<CrystalNeighborBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        self.values.push(value);
        if self.values.len() == self.capacity {
            self.flush(context, visitor)?;
        }
        Ok(())
    }

    fn flush(
        &mut self,
        context: &ExecutionContext,
        visitor: &mut impl FnMut(BatchLease<CrystalNeighborBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        let Some(reservation) = self.reservation.take() else {
            return Err(invariant());
        };
        let values = std::mem::take(&mut self.values);
        let lease = BatchLease::try_from_reservation(
            CrystalNeighborBatch { neighbors: values },
            reservation,
        )
        .map_err(|error| resource_error(&error.to_string()))?;
        visitor(lease)?;
        let bytes = self
            .capacity
            .checked_mul(std::mem::size_of::<CrystalNeighbor>())
            .ok_or_else(search_limit)?;
        self.reservation = Some(reserve(context, bytes)?);
        self.values = Vec::with_capacity(self.capacity);
        Ok(())
    }

    fn finish(
        mut self,
        visitor: &mut impl FnMut(BatchLease<CrystalNeighborBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        if self.values.is_empty() {
            return Ok(());
        }
        let Some(reservation) = self.reservation.take() else {
            return Err(invariant());
        };
        let lease = BatchLease::try_from_reservation(
            CrystalNeighborBatch {
                neighbors: self.values,
            },
            reservation,
        )
        .map_err(|error| resource_error(&error.to_string()))?;
        visitor(lease)
    }
}

fn reserve(context: &ExecutionContext, bytes: usize) -> Result<MemoryReservation, Diagnostic> {
    context
        .try_reserve(bytes)
        .map_err(|error| resource_error(&error.to_string()))
}

fn resource_error(reason: &str) -> Diagnostic {
    Diagnostic::new(Code::E1902).with_context("reason", reason)
}

fn search_limit() -> Diagnostic {
    Diagnostic::new(Code::E6017)
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

#[cfg(test)]
#[path = "crystal_tests.rs"]
mod tests;
