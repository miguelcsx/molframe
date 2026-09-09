//! Bounded crystal-image production without result-size materialisation.

use crate::SymmetrySet;
use crate::crystal_images::{CrystalImage, PreparedImages, prepare_images, visit_prepared_images};
use pdbiox_core::{
    Batch, BatchDemand, BatchLease, Code, Diagnostic, ExecutionContext, MemoryReservation,
    ModelIndex, Structure,
};

/// One caller-bounded run of explicit crystal images.
#[derive(Debug)]
pub struct CrystalImageBatch {
    images: Vec<CrystalImage>,
}

/// Resource and work controls for crystal-image batches.
#[derive(Clone, Copy, Debug)]
pub struct CrystalImageBatchOptions {
    /// Cartesian cutoff used to derive lattice images.
    pub cutoff: f64,
    /// Maximum candidate images examined.
    pub candidate_limit: usize,
    /// Maximum rows and retained bytes per emitted batch.
    pub demand: BatchDemand,
}

impl CrystalImageBatch {
    /// Images in deterministic generation order.
    #[must_use]
    pub fn images(&self) -> &[CrystalImage] {
        &self.images
    }
}

impl Batch for CrystalImageBatch {
    fn rows(&self) -> usize {
        self.images.len()
    }

    fn retained_bytes(&self) -> usize {
        self.images
            .capacity()
            .saturating_mul(std::mem::size_of::<CrystalImage>())
    }
}

/// Visits crystal images one at a time in stable generation order.
///
/// # Errors
///
/// Returns a crystal diagnostic for invalid input, cancellation, arithmetic
/// overflow, or work above `limit` candidate comparisons.
pub fn visit_crystal_images(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    limit: usize,
    context: &ExecutionContext,
    visitor: impl FnMut(CrystalImage) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let prepared = prepare_images(structure, symmetry, model, cutoff, context)?;
    visit_prepared_images(&prepared, symmetry, limit, context, visitor)
}

/// Visits budget-charged image batches in stable generation order.
///
/// The callback owns each lease for its duration. Keeping a lease alive through
/// a later pull consumes the same shared context capacity and therefore applies
/// backpressure to the complete execution graph.
///
/// # Errors
///
/// Returns a typed diagnostic when demand cannot hold one image, the context
/// cannot reserve the requested batch, or image generation fails.
pub fn crystal_image_batches(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    options: CrystalImageBatchOptions,
    context: &ExecutionContext,
    visitor: impl FnMut(BatchLease<CrystalImageBatch>) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let prepared = prepare_images(structure, symmetry, model, options.cutoff, context)?;
    visit_prepared_image_batches(
        &prepared,
        symmetry,
        options.candidate_limit,
        options.demand,
        context,
        visitor,
    )
}

pub(crate) fn visit_prepared_image_batches(
    prepared: &PreparedImages,
    symmetry: &SymmetrySet,
    limit: usize,
    demand: BatchDemand,
    context: &ExecutionContext,
    mut visitor: impl FnMut(BatchLease<CrystalImageBatch>) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let mut accumulator = ImageAccumulator::new(demand, context)?;
    visit_prepared_images(prepared, symmetry, limit, context, |image| {
        accumulator.push(image, context, &mut visitor)
    })?;
    accumulator.finish(&mut visitor)
}

struct ImageAccumulator {
    images: Vec<CrystalImage>,
    reservation: Option<MemoryReservation>,
    capacity: usize,
}

impl ImageAccumulator {
    fn new(demand: BatchDemand, context: &ExecutionContext) -> Result<Self, Diagnostic> {
        let width = std::mem::size_of::<CrystalImage>();
        let capacity = demand.max_rows.min(demand.max_bytes / width);
        if capacity == 0 {
            return Err(resource_error(
                "batch demand cannot retain one crystal image",
            ));
        }
        let bytes = capacity.checked_mul(width).ok_or_else(|| {
            resource_error("crystal image batch capacity exceeds addressable memory")
        })?;
        let reservation = reserve(context, bytes)?;
        Ok(Self {
            images: Vec::with_capacity(capacity),
            reservation: Some(reservation),
            capacity,
        })
    }

    fn push(
        &mut self,
        image: CrystalImage,
        context: &ExecutionContext,
        visitor: &mut impl FnMut(BatchLease<CrystalImageBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        self.images.push(image);
        if self.images.len() == self.capacity {
            self.flush(context, visitor)?;
        }
        Ok(())
    }

    fn flush(
        &mut self,
        context: &ExecutionContext,
        visitor: &mut impl FnMut(BatchLease<CrystalImageBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        if self.images.is_empty() {
            return Ok(());
        }
        let next_bytes = self
            .capacity
            .checked_mul(std::mem::size_of::<CrystalImage>())
            .ok_or_else(|| resource_error("crystal image batch capacity overflow"))?;
        let Some(reservation) = self.reservation.take() else {
            return Err(resource_error("crystal image reservation is absent"));
        };
        let images = std::mem::take(&mut self.images);
        let lease = BatchLease::try_from_reservation(CrystalImageBatch { images }, reservation)
            .map_err(|error| resource_error(&error.to_string()))?;
        visitor(lease)?;
        self.reservation = Some(reserve(context, next_bytes)?);
        self.images = Vec::with_capacity(self.capacity);
        Ok(())
    }

    fn finish(
        mut self,
        visitor: &mut impl FnMut(BatchLease<CrystalImageBatch>) -> Result<(), Diagnostic>,
    ) -> Result<(), Diagnostic> {
        if self.images.is_empty() {
            return Ok(());
        }
        let Some(reservation) = self.reservation.take() else {
            return Err(resource_error("crystal image reservation is absent"));
        };
        let lease = BatchLease::try_from_reservation(
            CrystalImageBatch {
                images: self.images,
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
