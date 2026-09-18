use super::*;
use molframe_core::{InputBuffer, ScratchPolicy};

const TWO_ATOMS: &str = concat!(
    "ATOM      1  N   GLY A   1      11.104  13.207  10.111  1.00 20.00           N  \n",
    "ATOM      2  CA  GLY A   1      12.104  13.207  10.111  1.00 21.00           C  \n",
    "END\n",
);

fn context() -> ExecutionContext {
    ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context")
}

#[test]
fn batches_advance_once_and_keep_global_rows_across_a_split_residue() {
    let source = InputBuffer::from_bytes(TWO_ATOMS.as_bytes().to_vec());
    let mut reader = PdbBatchSource::new(
        source,
        ReadOptions::new(),
        DatasetId::new(7),
        ChunkId::new(11),
        LogicalRow::new(u64::from(u32::MAX) + 4),
        4096,
    )
    .expect("valid source");
    let context = context();
    let demand = BatchDemand::new(1, 128 * 1024);

    let Backpressure::Ready(first) = reader.next_batch(demand, &context).expect("first batch")
    else {
        panic!("first batch should be ready");
    };
    assert_eq!(first.batch().descriptor().chunk(), ChunkId::new(11));
    assert_eq!(
        first.batch().descriptor().logical_start().get(),
        4_294_967_299
    );
    assert_eq!(first.batch().continuity().after, ContinuityLevel::Residue);
    drop(first);

    let Backpressure::Ready(second) = reader.next_batch(demand, &context).expect("second batch")
    else {
        panic!("second batch should be ready");
    };
    assert_eq!(second.batch().descriptor().chunk(), ChunkId::new(12));
    assert_eq!(
        second.batch().descriptor().logical_start().get(),
        4_294_967_300
    );
    assert_eq!(second.batch().continuity().before, ContinuityLevel::Residue);
    for (actual, expected) in second.batch().positions()[0]
        .iter()
        .zip([12.104, 13.207, 10.111])
    {
        assert!((*actual - expected).abs() <= f32::EPSILON);
    }
}

#[test]
fn zero_demand_does_not_advance_the_source() {
    let source = InputBuffer::from_bytes(TWO_ATOMS.as_bytes().to_vec());
    let mut reader = PdbBatchSource::new(
        source,
        ReadOptions::new(),
        DatasetId::new(1),
        ChunkId::new(0),
        LogicalRow::new(0),
        4096,
    )
    .expect("valid source");
    let context = context();
    assert!(matches!(
        reader
            .next_batch(BatchDemand::new(0, 0), &context)
            .expect("backpressure"),
        Backpressure::Pending
    ));
    assert!(matches!(
        reader
            .next_batch(BatchDemand::new(2, 128 * 1024), &context)
            .expect("unchanged source"),
        Backpressure::Ready(_)
    ));
}
