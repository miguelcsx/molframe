use super::*;
use crate::chunk::AtomRecord;
use crate::provider::{LocalRow, PropertyValue};
use crate::structure::CoordinateStore;
use crate::{
    AltId, AnnotationColumn, AtomAnnotation, AtomAnnotations, ChunkBuilder, Element,
    OptionalSymbol, Presence, ResidueIndex, StructureData, SymbolId,
};
use num_traits::ToPrimitive;

fn structure() -> Structure {
    let mut builder = ChunkBuilder::with_target(2);
    for atom in 0..5u32 {
        builder.push(AtomRecord {
            position: Some([atom.to_f32().expect("small atom ordinal"), 0.0, 0.0]),
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
    let (chunks, coordinates) = builder.finish();
    let mut annotations = AtomAnnotations::default();
    let property =
        AnnotationColumn::from_values(vec![0.0, 1.0, 2.0, 3.0, 4.0]).expect("small property");
    let _ = annotations.insert("score", AtomAnnotation::Real(property));
    let mut data = StructureData::empty();
    data.chunks = Arc::new(chunks);
    data.coords = CoordinateStore::Single(coordinates);
    data.annotations = annotations;
    Structure::new(data)
}

#[test]
fn structure_and_frame_providers_reuse_the_same_coordinate_storage() {
    let structure = structure();
    let pointer = structure.positions().as_ptr();
    let structures =
        StructureChunkProvider::new(DatasetId::new(7), ChunkId::new(100), structure.clone())
            .expect("structure provider");
    let frames = FrameChunkProvider::new(
        DatasetId::new(8),
        ChunkId::new(200),
        structure,
        crate::ModelIndex::new(0),
    )
    .expect("frame provider");

    let structure_chunk = structures.chunk(ChunkId::new(100)).expect("first chunk");
    let frame_chunk = frames.chunk(ChunkId::new(200)).expect("first frame chunk");
    assert_eq!(structure_chunk.positions().as_ptr(), pointer);
    assert_eq!(frame_chunk.positions().as_ptr(), pointer);
}

#[test]
fn property_chunks_read_native_columns_without_materialising_values() {
    let structure = structure();
    let source_pointer = match structure.annotations().get("score") {
        Some(AtomAnnotation::Real(column)) => column.values().as_ptr(),
        _ => panic!("real property expected"),
    };
    let provider =
        PropertyChunkProvider::new(DatasetId::new(9), ChunkId::new(300), structure, "score")
            .expect("property provider");
    let chunk = provider
        .chunk(ChunkId::new(300))
        .expect("first property chunk");
    let values = chunk.reals().expect("real values");

    assert_eq!(values.as_ptr(), source_pointer);
    assert_eq!(values, &[0.0, 1.0]);
    assert_eq!(
        chunk.value(LocalRow::new(1)).expect("local value"),
        PropertyValue::Real(1.0, Presence::Present)
    );
}
