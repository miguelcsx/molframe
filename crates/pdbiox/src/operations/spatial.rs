//! Typed spatial-search requests for the native facade plan.

use super::requests::{CoordinateInput, ExecutionPlanError, PlanInput};
use super::spatial_cache::SpatialContext;
use pdbiox_core::ExecutionContext;
use pdbiox_core::selection::AtomSelection;

/// A reusable fixed-radius spatial request over one borrowed coordinate array.
#[derive(Clone, Debug)]
pub enum SpatialRequest {
    /// Enumerates unique pairs within the cutoff.
    NeighborPairs {
        /// Coordinate-array slot.
        positions: usize,
        /// Left selection.
        left: AtomSelection,
        /// Right selection.
        right: AtomSelection,
        /// Distance cutoff in ångström.
        cutoff: f32,
        /// Complete backend planning profile.
        options: pdbiox_spatial::SpatialSearchOptions,
        /// Optional periodic unit-cell geometry.
        periodic: Option<pdbiox_spatial::PeriodicBox>,
    },
    /// Selects query atoms within the cutoff of a target selection.
    AtomsWithin {
        /// Coordinate-array slot.
        positions: usize,
        /// Query selection.
        query: AtomSelection,
        /// Target selection.
        target: AtomSelection,
        /// Distance cutoff in ångström.
        cutoff: f32,
        /// Complete backend planning profile.
        options: pdbiox_spatial::SpatialSearchOptions,
        /// Optional periodic unit-cell geometry.
        periodic: Option<pdbiox_spatial::PeriodicBox>,
    },
}

/// Results emitted by a spatial request.
#[derive(Clone, Debug)]
pub enum SpatialValue {
    /// Sorted fixed-radius pairs.
    NeighborPairs(Vec<pdbiox_spatial::NeighborPair>),
    /// Sorted atom indices selected by a within query.
    AtomsWithin(AtomSelection),
}

pub(super) fn execute(
    operation: &str,
    request: &SpatialRequest,
    input: PlanInput<'_>,
    context: Option<&mut SpatialContext<'_>>,
    execution: &ExecutionContext,
) -> Result<SpatialValue, ExecutionPlanError> {
    match request {
        SpatialRequest::NeighborPairs {
            positions,
            left,
            right,
            cutoff,
            options,
            periodic,
        } => {
            let result = match context {
                Some(context) => context.pairs(left, right, *cutoff, execution),
                None => pdbiox_spatial::pairs_within_with_options(
                    coordinates(operation, input.arrays, *positions)?,
                    left,
                    right,
                    *cutoff,
                    *options,
                    periodic.as_ref(),
                    execution,
                ),
            };
            result
                .map(SpatialValue::NeighborPairs)
                .map_err(|error| ExecutionPlanError::SpatialKernel(error.to_string().into()))
        }
        SpatialRequest::AtomsWithin {
            positions,
            query,
            target,
            cutoff,
            options,
            periodic,
        } => {
            let result = match context {
                Some(context) => context.within(query, target, *cutoff, execution),
                None => pdbiox_spatial::within_with_options(
                    coordinates(operation, input.arrays, *positions)?,
                    query,
                    target,
                    *cutoff,
                    *options,
                    periodic.as_ref(),
                    execution,
                ),
            };
            result
                .map(SpatialValue::AtomsWithin)
                .map_err(|error| ExecutionPlanError::SpatialKernel(error.to_string().into()))
        }
    }
}

fn coordinates<'a>(
    operation: &str,
    arrays: &'a [CoordinateInput<'a>],
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
