use super::*;

#[test]
fn edges_are_normalised_deduplicated_and_sorted() {
    let mut builder = BondTableBuilder::new();
    builder.push(BondRecord {
        atom_a: AtomIndex::new(2),
        atom_b: AtomIndex::new(0),
        order: BondOrder::Double,
        provenance: BondProvenance::File,
    });
    builder.push(BondRecord {
        atom_a: AtomIndex::new(0),
        atom_b: AtomIndex::new(2),
        order: BondOrder::Single,
        provenance: BondProvenance::InferredDistance,
    });
    builder.push(BondRecord {
        atom_a: AtomIndex::new(1),
        atom_b: AtomIndex::new(0),
        order: BondOrder::Single,
        provenance: BondProvenance::User,
    });
    let table = builder.finish();
    let records: Vec<_> = table.iter().collect();
    assert_eq!(records.len(), 2);
    assert_eq!((records[0].atom_a.get(), records[0].atom_b.get()), (0, 1));
    assert_eq!((records[1].atom_a.get(), records[1].atom_b.get()), (0, 2));
    assert_eq!(records[1].order, BondOrder::Double);
}

#[test]
fn adjacency_is_compact_bidirectional_and_stably_sorted() {
    let mut builder = BondTableBuilder::new();
    for endpoint in [3, 1, 2] {
        builder.push(BondRecord {
            atom_a: AtomIndex::new(0),
            atom_b: AtomIndex::new(endpoint),
            order: BondOrder::Single,
            provenance: BondProvenance::File,
        });
    }
    let table = builder.finish();
    assert_eq!(
        table.adjacency(4).neighbours(AtomIndex::new(0)),
        &[AtomIndex::new(1), AtomIndex::new(2), AtomIndex::new(3)]
    );
    assert_eq!(
        table.adjacency(4).neighbours(AtomIndex::new(2)),
        &[AtomIndex::new(0)]
    );
}
