use super::*;
use pdbiox::core::{Backpressure, BatchDemand, BatchSource};

#[test]
fn run_length_fixture_decodes_more_than_one_batch_exactly() {
    let rows = 200_000_u64;
    let context = ExecutionContext::default();
    let bytes = encoded_file(rows).expect("fixture");
    let mut source = pdbiox::bcif::BcifBatchSource::new(
        InputBuffer::from_bytes(bytes),
        ReadOptions::new(),
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
        64 * 1024,
        &context,
    )
    .expect("source");
    let mut actual = 0_u64;
    loop {
        match source
            .next_batch(BatchDemand::new(65_536, 16 * 1024 * 1024), &context)
            .expect("batch")
        {
            Backpressure::Ready(batch) => actual += batch.batch().models().len() as u64,
            Backpressure::Finished => break,
            Backpressure::Pending => panic!("unexpected backpressure"),
        }
    }
    assert_eq!(actual, rows);
}

#[test]
fn compact_fixture_represents_logical_rows_beyond_u32() {
    let rows = u64::from(u32::MAX) + 17;
    let bytes = encoded_file(rows).expect("fixture");
    assert!(bytes.len() < 4096);
}
