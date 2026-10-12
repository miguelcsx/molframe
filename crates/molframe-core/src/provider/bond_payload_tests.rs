use super::*;
use crate::StructureData;

#[test]
fn unavailable_topology_is_rejected_explicitly() {
    let descriptor = ChunkDescriptor::new(
        DatasetId::new(1),
        super::super::ChunkId::new(2),
        LogicalRow::new(0),
        1,
    )
    .expect("valid descriptor");
    let result = BondChunk::shared(
        descriptor,
        Structure::new(StructureData::empty()),
        0..1,
        DatasetId::new(3),
        LogicalRow::new(0),
    );

    assert!(matches!(
        result,
        Err(ProviderError::BondTopologyUnavailable)
    ));
}

/// The projected record carries the connection's own classification.
///
/// A bond read across a chunk boundary has to mean the same thing as the bond it
/// came from. Dropping the class here would make metal coordination and covalent
/// connectivity indistinguishable to every provider consumer, which is the whole
/// reason the class is stored separately from the order.
#[test]
fn a_projected_record_retains_the_connection_class() {
    use crate::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};

    // A class that is not canonical, so it occupies a dynamic dictionary slot and
    // a dropped field cannot pass by coincidence.
    let mut data = StructureData::empty();
    let Ok(kind) = data.dictionary.intern("source-class-alpha") else {
        panic!("interning failed")
    };
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: crate::AtomIndex::new(0),
        atom_b: crate::AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
        kind: crate::optional::OptionalSymbol::some(kind),
    });
    data.bonds = bonds.finish();
    let structure = Structure::new(data);

    let descriptor = ChunkDescriptor::new(
        DatasetId::new(1),
        super::super::ChunkId::new(2),
        LogicalRow::new(0),
        1,
    )
    .expect("valid descriptor");
    let chunk = BondChunk::shared(
        descriptor,
        structure.clone(),
        0..1,
        DatasetId::new(3),
        LogicalRow::new(0),
    )
    .expect("the bond table is available and the span is in range");

    let record = chunk
        .record(LocalRow::new(0))
        .expect("row zero is in range");
    assert_eq!(record.order, BondOrder::Single);
    assert_eq!(record.provenance, BondProvenance::File);
    assert_eq!(
        record.kind,
        crate::optional::OptionalSymbol::some(kind),
        "the connection class must survive the projection"
    );
    // And it still names the same string in the structure it came from.
    assert_eq!(
        chunk.structure().resolve(kind).map(str::to_owned),
        Some("source-class-alpha".to_owned())
    );
}
