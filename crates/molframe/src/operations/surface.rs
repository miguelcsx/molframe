//! Typed surface requests executed by the native facade plan.

use super::plan::inputs::PlanInput;
use super::plan::value::ExecutionPlanError;
use molframe_core::ExecutionContext;

/// A reusable surface calculation over borrowed plan inputs.
#[derive(Clone, Debug)]
pub enum SurfaceRequest {
    /// Per-atom Shrake–Rupley solvent-accessible areas.
    SolventAccessibleSurface {
        /// Coordinate-array slot.
        positions: usize,
        /// Per-atom radius-array slot.
        radii: usize,
        /// Probe radius in ångström.
        probe: f32,
        /// Number of deterministic sphere samples.
        sample_points: u16,
    },
    /// Total accessible areas and buried interface area for a two-way mask.
    BuriedSurface {
        /// Coordinate-array slot.
        positions: usize,
        /// Per-atom radius-array slot.
        radii: usize,
        /// Boolean mask-array slot; `true` selects the first group.
        first: usize,
        /// Probe radius in ångström.
        probe: f32,
        /// Number of deterministic sphere samples.
        sample_points: u16,
    },
}

/// Results emitted by a surface request.
#[derive(Clone, Debug)]
pub enum SurfaceValue {
    /// One accessible area per input atom.
    SolventAccessibleSurface(Vec<f64>),
    /// Accessible areas and the total buried area.
    BuriedSurface(molframe_surface::BuriedSurface),
}

pub(super) fn execute(
    operation: &str,
    request: &SurfaceRequest,
    input: PlanInput<'_>,
    context: &ExecutionContext,
) -> Result<SurfaceValue, ExecutionPlanError> {
    match request {
        SurfaceRequest::SolventAccessibleSurface {
            positions,
            radii,
            probe,
            sample_points,
        } => molframe_surface::shrake_rupley(
            coordinates(operation, input.arrays, *positions)?,
            floats(operation, input.floats, *radii)?,
            *probe,
            *sample_points,
            context,
        )
        .map(SurfaceValue::SolventAccessibleSurface)
        .map_err(|error| ExecutionPlanError::SurfaceKernel(error.to_string().into())),
        SurfaceRequest::BuriedSurface {
            positions,
            radii,
            first,
            probe,
            sample_points,
        } => molframe_surface::buried_surface(
            coordinates(operation, input.arrays, *positions)?,
            floats(operation, input.floats, *radii)?,
            *probe,
            *sample_points,
            masks(operation, input.masks, *first)?,
            context,
        )
        .map(SurfaceValue::BuriedSurface)
        .map_err(|error| ExecutionPlanError::SurfaceKernel(error.to_string().into())),
    }
}

fn coordinates<'a>(
    operation: &str,
    arrays: &'a [super::plan::inputs::CoordinateInput<'a>],
    slot: usize,
) -> Result<&'a [[f32; 3]], ExecutionPlanError> {
    arrays
        .get(slot)
        .map(|input| input.positions)
        .ok_or_else(|| ExecutionPlanError::ArraySlot {
            operation: operation.into(),
            slot,
        })
}

fn floats<'a>(
    operation: &str,
    values: &'a [super::plan::inputs::FloatInput<'a>],
    slot: usize,
) -> Result<&'a [f32], ExecutionPlanError> {
    values
        .get(slot)
        .map(|input| input.values)
        .ok_or_else(|| ExecutionPlanError::FloatSlot {
            operation: operation.into(),
            slot,
        })
}

fn masks<'a>(
    operation: &str,
    values: &'a [super::plan::inputs::MaskInput<'a>],
    slot: usize,
) -> Result<&'a [bool], ExecutionPlanError> {
    values
        .get(slot)
        .map(|input| input.values)
        .ok_or_else(|| ExecutionPlanError::MaskSlot {
            operation: operation.into(),
            slot,
        })
}
