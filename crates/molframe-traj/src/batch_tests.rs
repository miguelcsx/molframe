use super::*;
use crate::{Frame, StreamingReader};
use molframe_core::{MemoryBudget, ScratchPolicy};

fn context(bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(bytes).expect("non-zero budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context")
}

#[test]
fn pull_honours_demand_and_stable_identity() {
    let reader = StreamingReader::new(
        vec![Frame {
            positions: vec![[1.0, 2.0, 3.0]],
        }]
        .into_iter(),
        1,
    );
    let mut source = TrajectoryBatchSource::new(
        reader,
        DatasetId::new(7),
        ChunkId::new(9),
        LogicalRow::new(u64::from(u32::MAX) + 1),
    );
    let execution = context(4096);
    assert!(matches!(
        source
            .next_batch(BatchDemand::new(0, 4096), &execution)
            .expect("pending pull"),
        Backpressure::Pending
    ));
    let ready = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("frame pull");
    let Backpressure::Ready(lease) = ready else {
        panic!("expected one frame");
    };
    assert_eq!(lease.batch().dataset(), DatasetId::new(7));
    assert_eq!(lease.batch().chunk(), ChunkId::new(9));
    assert_eq!(
        lease.batch().logical_row(),
        LogicalRow::new(u64::from(u32::MAX) + 1)
    );
    let timestep = lease.batch().timestep().expect("live timestep");
    assert_eq!(timestep.positions, [[1.0, 2.0, 3.0]]);
}

#[test]
fn insufficient_context_capacity_does_not_consume_input() {
    let reader = StreamingReader::new(
        vec![Frame {
            positions: vec![[1.0, 2.0, 3.0]],
        }]
        .into_iter(),
        1,
    );
    let mut source = TrajectoryBatchSource::new(
        reader,
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
    );
    let execution = context(128);
    assert!(matches!(
        source
            .next_batch(BatchDemand::new(1, 256), &execution)
            .expect("pending pull"),
        Backpressure::Pending
    ));
}

#[test]
fn cancellation_stops_before_advancing_the_reader() {
    let reader = StreamingReader::new(
        vec![Frame {
            positions: vec![[1.0, 2.0, 3.0]],
        }]
        .into_iter(),
        1,
    );
    let mut source = TrajectoryBatchSource::new(
        reader,
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
    );
    let execution = context(4096);
    execution.cancellation().cancel();
    assert_eq!(
        source
            .next_batch(BatchDemand::new(1, 4096), &execution)
            .expect_err("cancelled pull"),
        TrajectoryError::Cancelled
    );
}

#[test]
fn recycled_frame_capacity_stays_charged_and_is_reused() {
    let execution = context(4096);
    let reader = StreamingReader::new(
        (0..2).map(|_| Frame {
            positions: vec![[0.0; 3]; 16],
        }),
        16,
    );
    let mut source = TrajectoryBatchSource::new(
        reader,
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
    );
    let Backpressure::Ready(first) = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("first")
    else {
        panic!("ready")
    };
    let pointer = first.batch().timestep().expect("frame").positions.as_ptr();
    let charged = execution.reserved_bytes();
    drop(first);
    assert_eq!(execution.reserved_bytes(), charged);
    assert_eq!(execution.live_batches(), 0);
    let Backpressure::Ready(second) = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("second")
    else {
        panic!("ready")
    };
    assert_eq!(
        second.batch().timestep().expect("frame").positions.as_ptr(),
        pointer
    );
    drop(second);
    drop(source);
    assert_eq!(execution.reserved_bytes(), 0);
}

#[test]
fn recycling_retains_only_one_idle_frame() {
    let execution = context(8192);
    let reader = StreamingReader::new(
        (0..2).map(|_| Frame {
            positions: vec![[0.0; 3]; 16],
        }),
        16,
    );
    let mut source = TrajectoryBatchSource::new(
        reader,
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
    );
    let Backpressure::Ready(first) = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("first")
    else {
        panic!("ready")
    };
    let charged = execution.reserved_bytes();
    let second = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("second");
    assert_eq!(execution.reserved_bytes(), charged * 2);
    drop(first);
    drop(second);
    assert_eq!(execution.reserved_bytes(), charged);
    drop(source);
    assert_eq!(execution.reserved_bytes(), 0);
}

#[test]
fn exported_lease_metadata_is_reserved_before_input_is_consumed() {
    let mut source = TrajectoryBatchSource::new(
        StreamingReader::new(
            std::iter::once(Frame {
                positions: vec![[1.0; 3]],
            }),
            1,
        ),
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
    );
    let execution = context(4096);
    let only_payload = size_of::<Timestep>() + size_of::<[f32; 3]>();
    assert!(matches!(
        source.next_batch(BatchDemand::new(1, only_payload), &execution),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
    assert_eq!(execution.reserved_bytes(), 0);
    let Backpressure::Ready(lease) = source
        .next_batch(BatchDemand::new(1, 4096), &execution)
        .expect("retry without advancing")
    else {
        panic!("frame expected");
    };
    assert_eq!(lease.batch().timestep().expect("frame").frame, 0);
    assert!(lease.batch().retained_bytes() >= only_payload + LEASE_OVERHEAD);
}
