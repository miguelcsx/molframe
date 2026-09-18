use super::*;
use crate::{
    ChunkId, DatasetId, LogicalRow, ScratchPolicy, StructureAtomRecord, StructureBatchBuilder,
};

#[derive(Default)]
struct TwoBatches {
    next: usize,
}

impl BatchSource for TwoBatches {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<crate::BatchLease<Self::Batch>>, Self::Error> {
        if !demand.can_accept_work() {
            return Ok(Backpressure::Pending);
        }
        if self.next == 2 {
            return Ok(Backpressure::Finished);
        }
        let mut builder = StructureBatchBuilder::new(
            DatasetId::new(0),
            ChunkId::new(self.next as u64),
            LogicalRow::new(self.next as u64),
            1,
        )?;
        let coordinate = [0.0_f32, 1.0][self.next];
        let atom_site_id = [1_u32, 2][self.next];
        builder.push(StructureAtomRecord {
            model: 1,
            chain: "A",
            component: "GLY",
            sequence: 1,
            insertion: "",
            atom: if self.next == 0 { "N" } else { "CA" },
            alternate: "",
            element: Element::CARBON,
            position: Some([coordinate, 0.0, 0.0]),
            occupancy: (1.0, crate::Presence::Present),
            b_factor: (0.0, crate::Presence::Present),
            formal_charge: (0, crate::Presence::Inapplicable),
            atom_site_id,
            heterogen: false,
        })?;
        self.next += 1;
        Ok(Backpressure::Ready(crate::BatchLease::try_new(
            builder.finish()?.into_owned(),
            context,
        )?))
    }
}

#[test]
fn collector_appends_batches_without_structure_merge() {
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let (structure, findings) = collect_structure(&mut TwoBatches::default(), &context)
        .expect("bounded collection succeeds");
    assert_eq!(structure.atom_count(), 2);
    assert_eq!(structure.residue_count(), 1);
    assert!(findings.is_empty());
}

#[test]
fn collector_rejects_growth_before_a_tiny_budget_is_exceeded() {
    let budget = crate::MemoryBudget::new(64 * 1024).expect("non-zero budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let error = collect_structure(&mut TwoBatches::default(), &context)
        .expect_err("collector base consumes the complete budget");
    assert!(matches!(
        error,
        StructureBatchError::RecordExceedsBudget { .. }
    ));
}

#[test]
fn collected_snapshot_and_shared_aliases_retain_the_reservation() {
    let context = ExecutionContext::default();
    let (structure, _) = collect_structure(&mut TwoBatches::default(), &context).expect("collect");
    let charge = context.reserved_bytes();
    assert!(charge >= COLLECTOR_BASE_BYTES + 2 * COLLECTOR_BYTES_PER_ROW);
    let alias = structure.clone();
    let shared = structure.shared();
    drop(structure);
    drop(alias);
    assert_eq!(context.reserved_bytes(), charge);
    drop(shared);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn an_independent_ragged_model_keeps_the_shared_allocation_charge() {
    let context = ExecutionContext::default();
    let reservation = std::sync::Arc::new(context.try_reserve(1024).expect("reserve"));
    let mut data = StructureData::empty();
    data.coords = CoordinateStore::Ragged {
        models: vec![
            Structure::new(StructureData::empty()),
            Structure::new(StructureData::empty()),
        ],
    };
    let structure = Structure::new(data).retain_reservation(reservation);
    let model = structure.ragged_models().expect("models")[0].clone();
    drop(structure);
    assert_eq!(context.reserved_bytes(), 1024);
    drop(model);
    assert_eq!(context.reserved_bytes(), 0);
}
