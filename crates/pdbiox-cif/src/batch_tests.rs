use super::*;
use pdbiox_core::{ByteWindow, Diagnostic, InputBuffer, ScratchPolicy};
use std::cell::Cell;
use std::rc::Rc;

const TWO_ATOMS: &str = "data_x\nloop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.pdbx_formal_charge\n\
_atom_site.pdbx_PDB_model_num\n\
ATOM 1 N N . GLY A 1 ? 1.0 2.0 3.0 1.0 10.0 ? 1\n\
ATOM 2 C CA . GLY A 1 ? 4.0 5.0 6.0 1.0 11.0 ? 1\n";

fn context() -> ExecutionContext {
    ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context")
}

struct CountingSource {
    input: InputBuffer,
    calls: Rc<Cell<usize>>,
}

impl SourceBytes for CountingSource {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic> {
        self.calls.set(self.calls.get().saturating_add(1));
        self.input.window(start, len)
    }

    fn len_hint(&self) -> Option<u64> {
        self.input.len_hint()
    }
}

#[test]
fn cif_rows_cross_batch_and_u32_boundaries_without_losing_identity() {
    let context = context();
    let source = InputBuffer::from_bytes(TWO_ATOMS.as_bytes().to_vec());
    let mut reader = MmcifBatchSource::new(
        source,
        ReadOptions::new(),
        DatasetId::new(3),
        ChunkId::new(9),
        LogicalRow::new(u64::from(u32::MAX)),
        128,
        &context,
    )
    .expect("valid source");
    let demand = BatchDemand::new(1, 128 * 1024);

    let Backpressure::Ready(first) = reader.next_batch(demand, &context).expect("first batch")
    else {
        panic!("first batch should be ready");
    };
    assert_eq!(
        first.batch().descriptor().logical_start().get(),
        4_294_967_295
    );
    assert!(
        first.batch().positions()[0]
            .iter()
            .zip([1.0, 2.0, 3.0])
            .all(|(actual, expected)| (*actual - expected).abs() <= f32::EPSILON)
    );
    drop(first);

    let Backpressure::Ready(second) = reader.next_batch(demand, &context).expect("second batch")
    else {
        panic!("second batch should be ready");
    };
    assert_eq!(
        second.batch().descriptor().logical_start().get(),
        4_294_967_296
    );
    assert!(
        second.batch().positions()[0]
            .iter()
            .zip([4.0, 5.0, 6.0])
            .all(|(actual, expected)| (*actual - expected).abs() <= f32::EPSILON)
    );
}

#[test]
fn a_comment_longer_than_one_window_is_discarded_without_becoming_a_token() {
    let mut text = String::from("#");
    text.push_str(&"x".repeat(300));
    text.push('\n');
    text.push_str(TWO_ATOMS);
    let context = context();
    let source = InputBuffer::from_bytes(text.into_bytes());
    let mut reader = MmcifBatchSource::new(
        source,
        ReadOptions::new(),
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
        128,
        &context,
    )
    .expect("valid source");
    let Backpressure::Ready(batch) = reader
        .next_batch(BatchDemand::new(2, 128 * 1024), &context)
        .expect("batch after long comment")
    else {
        panic!("coordinate batch should be ready");
    };
    assert_eq!(batch.batch().rows(), 2);
}

#[test]
fn tokenization_reuses_each_source_window_instead_of_rereading_it_per_token() {
    let calls = Rc::new(Cell::new(0));
    let source = CountingSource {
        input: InputBuffer::from_bytes(TWO_ATOMS.as_bytes().to_vec()),
        calls: Rc::clone(&calls),
    };
    let context = context();
    let mut reader = MmcifBatchSource::new(
        source,
        ReadOptions::new(),
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
        128,
        &context,
    )
    .expect("valid source");
    let Backpressure::Ready(batch) = reader
        .next_batch(BatchDemand::new(2, 128 * 1024), &context)
        .expect("coordinate batch")
    else {
        panic!("coordinate batch should be ready");
    };
    assert_eq!(batch.batch().rows(), 2);
    assert!(calls.get() <= 8, "source windows: {}", calls.get());
}

#[test]
fn demand_rollback_restarts_from_the_same_complete_row() {
    let long_atom = "A".repeat(80_000);
    let text = TWO_ATOMS.replacen("ATOM 1 N N .", &format!("ATOM 1 N {long_atom} ."), 1);
    let context = context();
    let mut reader = MmcifBatchSource::new(
        InputBuffer::from_bytes(text.into_bytes()),
        ReadOptions::new(),
        DatasetId::new(0),
        ChunkId::new(0),
        LogicalRow::new(0),
        100_000,
        &context,
    )
    .expect("valid source");
    let error = reader
        .next_batch(BatchDemand::new(1, 66_000), &context)
        .expect_err("long identifier exceeds the first demand");
    assert!(matches!(error, StructureBatchError::DemandTooSmall { .. }));
    let Backpressure::Ready(batch) = reader
        .next_batch(BatchDemand::new(1, 256 * 1024), &context)
        .expect("retry")
    else {
        panic!("retry should return the first row");
    };
    assert_eq!(
        batch.batch().descriptor().logical_start(),
        LogicalRow::new(0)
    );
    let atom = batch.batch().dictionary().resolve(batch.batch().atoms()[0]);
    assert_eq!(atom, Some(long_atom.as_str()));
}
