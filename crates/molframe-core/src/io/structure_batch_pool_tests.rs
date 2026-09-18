use super::*;
use crate::{
    ChunkId, DatasetId, Element, LogicalRow, Presence, ScratchPolicy, StructureAtomRecord,
};

#[test]
fn a_released_lease_reuses_columns_without_hiding_their_reservation() {
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let pool = StructureBatchPool::new();
    let descriptor =
        ChunkDescriptor::new(DatasetId::new(1), ChunkId::new(0), LogicalRow::new(0), 4)
            .expect("descriptor");
    let (mut builder, reservation) = pool
        .acquire(descriptor, 4096, 4096, &context)
        .expect("acquire")
        .expect("capacity");
    builder.push(atom("A")).expect("first row");
    let batch = builder.finish().expect("first batch");
    let pointer = batch.batch().models().as_ptr();
    let lease = pool.lease(batch, reservation).expect("first lease");
    drop(lease);
    assert!(context.reserved_bytes() > 0);

    let descriptor =
        ChunkDescriptor::new(DatasetId::new(1), ChunkId::new(1), LogicalRow::new(1), 4)
            .expect("descriptor");
    let (mut builder, reservation) = pool
        .acquire(descriptor, 1, 4096, &context)
        .expect("reacquire")
        .expect("capacity");
    builder.push(atom("B")).expect("second row");
    let batch = builder.finish().expect("second batch");
    assert_eq!(batch.batch().models().as_ptr(), pointer);
    assert_eq!(
        batch.batch().descriptor().logical_start(),
        LogicalRow::new(1)
    );
    let lease = pool.lease(batch, reservation).expect("second lease");
    drop(lease);
    drop(pool);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn an_undercharged_rejected_buffer_is_dropped_instead_of_hidden() {
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let pool = StructureBatchPool::new();
    let mut builder =
        StructureBatchBuilder::new(DatasetId::new(1), ChunkId::new(0), LogicalRow::new(0), 1)
            .expect("builder");
    builder.push(atom("A")).expect("row");
    let buffer = builder.finish().expect("batch");
    let reservation = context.try_reserve(1).expect("tiny reservation");
    pool.recycle(buffer, reservation);
    assert_eq!(context.reserved_bytes(), 0);
}

fn atom(chain: &str) -> StructureAtomRecord<'_> {
    StructureAtomRecord {
        model: 1,
        chain,
        component: "ALA",
        sequence: 1,
        insertion: "",
        atom: "CA",
        alternate: "",
        element: Element::CARBON,
        position: Some([1.0, 2.0, 3.0]),
        occupancy: (1.0, Presence::Present),
        b_factor: (10.0, Presence::Present),
        formal_charge: (0, Presence::Unknown),
        atom_site_id: 1,
        heterogen: false,
    }
}
