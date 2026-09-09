use super::*;
use crate::LocalRow;
use crate::chunk::AtomRecord;
use crate::structure::CoordinateStore;
use crate::{
    AltId, AtomIndex, BondOrder, BondProvenance, BondRecord, BondTableBuilder, ChunkBuilder,
    Element, OptionalSymbol, Presence, ResidueIndex, StructureData, SymbolId,
};
use std::sync::Arc;

fn bonded_structure() -> Structure {
    let mut atoms = ChunkBuilder::with_target(2);
    for atom in 0..4u32 {
        atoms.push(AtomRecord {
            position: Some([0.0, 0.0, 0.0]),
            element: Element::CARBON,
            atom_name: SymbolId::from_raw(0),
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue: ResidueIndex::new(atom),
            occupancy: (1.0, Presence::Present),
            b_factor: (10.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: atom,
        });
    }
    let (chunks, coordinates) = atoms.finish();
    let mut bonds = BondTableBuilder::new();
    for (atom_a, atom_b) in [(0, 1), (1, 2), (2, 3)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(atom_a),
            atom_b: AtomIndex::new(atom_b),
            order: BondOrder::Single,
            provenance: BondProvenance::File,
        });
    }
    let mut data = StructureData::empty();
    data.chunks = Arc::new(chunks);
    data.coords = CoordinateStore::Single(coordinates);
    data.bonds = bonds.finish();
    Structure::new(data)
}

#[test]
fn requested_chunks_project_global_endpoints_without_copying_coordinates() {
    let structure = bonded_structure();
    let coordinates = structure.positions().as_ptr();
    let atom_dataset = DatasetId::new(u64::from(u32::MAX) + 90);
    let atom_start = LogicalRow::new(u64::from(u32::MAX) + 500);
    let provider = BondChunkProvider::with_rows_per_chunk(
        DatasetId::new(70),
        ChunkId::new(u64::from(u32::MAX) + 10),
        atom_dataset,
        atom_start,
        structure,
        2,
    )
    .expect("valid bond provider");

    assert_eq!(provider.dataset().chunk_count(), 2);
    let chunk = provider
        .chunk(ChunkId::new(u64::from(u32::MAX) + 10))
        .expect("first bond chunk");
    let cross_chunk = chunk.record(LocalRow::new(1)).expect("cross-chunk bond");

    assert_eq!(cross_chunk.atom_a.dataset(), atom_dataset);
    assert_eq!(cross_chunk.atom_a.row().get(), atom_start.get() + 1);
    assert_eq!(cross_chunk.atom_b.row().get(), atom_start.get() + 2);
    assert_eq!(chunk.structure().positions().as_ptr(), coordinates);
}

#[test]
fn each_request_is_bounded_to_its_regular_chunk_extent() {
    let provider = BondChunkProvider::with_rows_per_chunk(
        DatasetId::new(80),
        ChunkId::new(900),
        DatasetId::new(81),
        LogicalRow::new(0),
        bonded_structure(),
        2,
    )
    .expect("valid bond provider");

    let first = provider.chunk(ChunkId::new(900)).expect("first chunk");
    let tail = provider.chunk(ChunkId::new(901)).expect("tail chunk");
    assert_eq!(first.descriptor().rows(), 2);
    assert_eq!(tail.descriptor().rows(), 1);
    assert!(matches!(
        provider.chunk(ChunkId::new(902)),
        Err(ProviderError::UnknownChunk { .. })
    ));
}

#[test]
fn global_atom_endpoint_overflow_is_typed_before_serving_chunks() {
    let result = BondChunkProvider::new(
        DatasetId::new(90),
        ChunkId::new(1),
        DatasetId::new(91),
        LogicalRow::new(u64::MAX - 1),
        bonded_structure(),
    );

    assert!(matches!(
        result,
        Err(ProviderError::IdentityOverflow {
            field: "atom endpoint range",
            ..
        })
    ));
}
