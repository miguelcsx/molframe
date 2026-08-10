//! Typed frame analyses with worker-count-independent reduction order.

use crate::{Timestep, Trajectory, TrajectoryError};
use pdbiox_core::contract::Analysis;
use std::ops::Range;

/// Fixed leaf size keeps floating-point grouping independent of worker count.
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
    workers: usize,
) -> Result<Analysis<A::Output>, A::Error> {
    if workers == 0 {
        return Err(TrajectoryError::InvalidWorkerCount.into());
    }
    if workers > 1 && !A::PARALLELIZABLE {
        return Err(TrajectoryError::AnalysisNotParallel.into());
    }
    let blocks = blocks(trajectory.len());
    let partials = if blocks.is_empty() {
        vec![analysis.prepare(trajectory)?]
    } else if workers == 1 {
        sequential_blocks(trajectory, analysis, &blocks)?
    } else {
        parallel_blocks(trajectory, analysis, &blocks, workers)?
    };
    analysis.conclude(analysis.merge(partials)?)
}

fn blocks(frame_count: usize) -> Vec<Range<usize>> {
    (0..frame_count)
        .step_by(CANONICAL_BLOCK_SIZE)
        .map(|start| start..(start + CANONICAL_BLOCK_SIZE).min(frame_count))
        .collect()
}

fn sequential_blocks<A: FrameAnalysis>(
    trajectory: &Trajectory,
    analysis: &A,
    blocks: &[Range<usize>],
) -> Result<Vec<A::Partial>, A::Error> {
    blocks
        .iter()
        .map(|range| process_block(trajectory, analysis, range.clone()))
        .collect()
}

fn parallel_blocks<A: FrameAnalysis>(
    trajectory: &Trajectory,
    analysis: &A,
    blocks: &[Range<usize>],
    workers: usize,
) -> Result<Vec<A::Partial>, A::Error> {
    let worker_count = workers.min(blocks.len());
    let gathered = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for worker in 0..worker_count {
            handles.push(scope.spawn(move || {
                let mut output = Vec::new();
                for block_index in (worker..blocks.len()).step_by(worker_count) {
                    output.push((
                        block_index,
                        process_block(trajectory, analysis, blocks[block_index].clone())?,
                    ));
                }
                Ok::<_, A::Error>(output)
            }));
        }
        let mut output = Vec::with_capacity(blocks.len());
        for handle in handles {
            output.extend(
                handle
                    .join()
                    .map_err(|_| A::Error::from(TrajectoryError::WorkerPanicked))??,
            );
        }
        Ok::<_, A::Error>(output)
    })?;
    let mut ordered = gathered;
    ordered.sort_unstable_by_key(|(block_index, _)| *block_index);
    Ok(ordered.into_iter().map(|(_, partial)| partial).collect())
}

fn process_block<A: FrameAnalysis>(
    trajectory: &Trajectory,
    analysis: &A,
    range: Range<usize>,
) -> Result<A::Partial, A::Error> {
    let mut partial = analysis.prepare(trajectory)?;
    let mut timestep = Timestep::default();
    for frame_index in range {
        let Some(frame) = trajectory.frame(frame_index) else {
            return Err(TrajectoryError::MissingFrame { index: frame_index }.into());
        };
        timestep.frame = frame_index;
        timestep.positions.clear();
        timestep.positions.extend_from_slice(&frame.positions);
        analysis.single_frame(&timestep, &mut partial)?;
    }
    Ok(partial)
}

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
