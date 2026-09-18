//! Typed frame analyses with worker-count-independent reduction order.

use crate::{Timestep, Trajectory, TrajectoryError};
use molframe_core::ExecutionContext;
use molframe_core::contract::Analysis;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use std::ops::Range;

/// Fixed leaf size keeps floating-point grouping independent of worker count.
///
/// Changing it changes the numerical result, which is why it is a named
/// constant rather than a tuning parameter.
const CANONICAL_BLOCK_SIZE: usize = 64;

/// Analysis lifecycle over trajectory frames.
pub trait FrameAnalysis: Sync {
    /// Final public value.
    type Output;
    /// Thread-local accumulation state.
    type Partial: Send;
    /// Error returned by the lifecycle, including trajectory execution errors.
    type Error: From<TrajectoryError> + Send;

    /// Stable analysis name recorded by higher-level provenance adapters.
    const NAME: &'static str;

    /// Whether independent canonical blocks may run concurrently.
    const PARALLELIZABLE: bool = false;

    /// Creates one empty block accumulator.
    ///
    /// # Errors
    ///
    /// Returns an analysis-specific trajectory error before frame processing.
    fn prepare(&self, trajectory: &Trajectory) -> Result<Self::Partial, Self::Error>;

    /// Accumulates one frame.
    ///
    /// # Errors
    ///
    /// Returns an analysis-specific trajectory error without producing a final result.
    fn single_frame(
        &self,
        timestep: &Timestep,
        partial: &mut Self::Partial,
        context: &ExecutionContext,
    ) -> Result<(), Self::Error>;

    /// Merges canonical blocks in ascending block order.
    ///
    /// # Errors
    ///
    /// Returns an analysis-specific error instead of a partial merge.
    fn merge(&self, partials: Vec<Self::Partial>) -> Result<Self::Partial, Self::Error>;

    /// Converts the merged accumulator to the governed result.
    ///
    /// # Errors
    ///
    /// Returns an analysis-specific error instead of a partial result.
    fn conclude(&self, partial: Self::Partial) -> Result<Analysis<Self::Output>, Self::Error>;
}

/// Executes fixed frame blocks and merges them in deterministic index order.
///
/// Worker count changes only which thread computes a block. Block boundaries
/// and the final merge order remain identical for 1, 2, 4 or more workers.
///
/// # Errors
///
/// Refuses zero workers, parallel execution without opt-in, worker panics and
/// any error returned by the analysis lifecycle.
pub fn run_analysis<A: FrameAnalysis>(
    trajectory: &Trajectory,
    analysis: &A,
    context: &ExecutionContext,
) -> Result<Analysis<A::Output>, A::Error> {
    if context.worker_budget() > 1 && !A::PARALLELIZABLE {
        return Err(TrajectoryError::AnalysisNotParallel.into());
    }
    let plan = BlockPlan::new(trajectory.len(), CANONICAL_BLOCK_SIZE);
    let partials = if plan.is_empty() {
        vec![analysis.prepare(trajectory)?]
    } else {
        let produced = map_blocks_in(plan, context, |_, range| {
            process_block(trajectory, analysis, range, context)
        })
        .map_err(|_| TrajectoryError::WorkerPanicked)?;
        produced.into_iter().collect::<Result<Vec<_>, A::Error>>()?
    };
    analysis.conclude(analysis.merge(partials)?)
}

fn process_block<A: FrameAnalysis>(
    trajectory: &Trajectory,
    analysis: &A,
    range: Range<usize>,
    context: &ExecutionContext,
) -> Result<A::Partial, A::Error> {
    let mut partial = analysis.prepare(trajectory)?;
    let mut timestep = Timestep::default();
    for frame_index in range {
        let Some(frame) = trajectory.frame(frame_index) else {
            return Err(TrajectoryError::MissingFrame { index: frame_index }.into());
        };
        timestep.frame = frame_index;
        timestep.positions.clear();
        timestep.positions.extend_from_slice(frame.positions);
        analysis.single_frame(&timestep, &mut partial, context)?;
    }
    Ok(partial)
}

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
