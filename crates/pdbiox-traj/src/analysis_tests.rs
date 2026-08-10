use super::{FrameAnalysis, run_analysis};
use crate::{Frame, Timestep, Trajectory, TrajectoryError};
use pdbiox_core::contract::{Analysis, AnalysisPolicy, Coverage};

struct SumX {
    policy: AnalysisPolicy,
}

impl FrameAnalysis for SumX {
    type Output = f64;
    type Partial = f64;
    type Error = TrajectoryError;

    const NAME: &'static str = "sum-x";

    const PARALLELIZABLE: bool = true;

    fn prepare(&self, _trajectory: &Trajectory) -> Result<Self::Partial, TrajectoryError> {
        Ok(0.0)
    }

    fn single_frame(
        &self,
        timestep: &Timestep,
        partial: &mut Self::Partial,
    ) -> Result<(), TrajectoryError> {
        *partial += f64::from(timestep.positions[0][0]);
        Ok(())
    }

    fn merge(&self, partials: Vec<Self::Partial>) -> Result<Self::Partial, TrajectoryError> {
        Ok(partials.into_iter().sum())
    }

    fn conclude(&self, partial: Self::Partial) -> Result<Analysis<Self::Output>, TrajectoryError> {
        Ok(Analysis::complete(
            partial,
            Coverage::complete(1),
            &self.policy,
        ))
    }
}

#[test]
fn merge_is_bit_identical_at_every_worker_count() {
    let frames = (0..1_025)
        .map(|index| Frame {
            positions: vec![[if index % 3 == 0 { 1.0e10 } else { 0.1 }, 0.0, 0.0]],
        })
        .collect();
    let trajectory = Trajectory::from_frames(frames);
    let analysis = SumX {
        policy: AnalysisPolicy::default(),
    };
    let outputs: Vec<_> = [1, 2, 4, 16]
        .into_iter()
        .map(|workers| {
            run_analysis(&trajectory, &analysis, workers).map_or_else(
                |error| panic!("analysis failed: {error}"),
                |result| result.value.to_bits(),
            )
        })
        .collect();
    assert!(outputs.windows(2).all(|pair| pair[0] == pair[1]));
}

struct Ordered;

impl FrameAnalysis for Ordered {
    type Output = ();
    type Partial = ();
    type Error = TrajectoryError;

    const NAME: &'static str = "ordered";

    fn prepare(&self, _trajectory: &Trajectory) -> Result<Self::Partial, TrajectoryError> {
        Ok(())
    }

    fn single_frame(
        &self,
        _timestep: &Timestep,
        _partial: &mut Self::Partial,
    ) -> Result<(), TrajectoryError> {
        Ok(())
    }

    fn merge(&self, _partials: Vec<Self::Partial>) -> Result<Self::Partial, TrajectoryError> {
        Ok(())
    }

    fn conclude(&self, partial: Self::Partial) -> Result<Analysis<Self::Output>, TrajectoryError> {
        Ok(Analysis::complete(
            partial,
            Coverage::complete(0),
            &AnalysisPolicy::default(),
        ))
    }
}

#[test]
fn order_dependent_analysis_cannot_accidentally_run_in_parallel() {
    let trajectory = Trajectory::default();
    assert!(matches!(
        run_analysis(&trajectory, &Ordered, 2),
        Err(TrajectoryError::AnalysisNotParallel)
    ));
}
