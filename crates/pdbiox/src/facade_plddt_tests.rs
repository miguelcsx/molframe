use super::*;

#[cfg(feature = "modelcif")]
#[test]
fn modelcif_plddt_becomes_a_shared_atom_property_without_changing_b_factors() {
    let source = include_str!("../../pdbiox-modelcif/tests/fixtures/modelcif_plddt.cif");
    let (structure, findings) = read_bytes(
        source.as_bytes().to_vec(),
        Some("modelcif_plddt.cif"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|values| panic!("read failed: {values:?}"));
    assert!(findings.is_empty());
    let b_factor = structure
        .atom(crate::AtomIndex::new(0))
        .and_then(crate::AtomRef::b_factor)
        .expect("fixture B-factor");
    assert!((b_factor - 11.0).abs() < f32::EPSILON);

    let source_pointer = match structure.annotations().get(crate::PLDDT_ANNOTATION) {
        Some(crate::AtomAnnotation::Real(column)) => column.values().as_ptr(),
        _ => panic!("lowered real pLDDT expected"),
    };
    let dataset = crate::DatasetId::new(u64::from(u32::MAX) + 9);
    let first_chunk = crate::ChunkId::new(u64::from(u32::MAX) + 17);
    let provider = crate::PropertyChunkProvider::plddt(dataset, first_chunk, structure)
        .expect("pLDDT provider");
    let chunk = provider.chunk(first_chunk).expect("first pLDDT chunk");
    let values = chunk.reals().expect("real pLDDT values");
    assert_eq!(values.as_ptr(), source_pointer);
    assert_eq!(chunk.descriptor().dataset(), dataset);
    assert_eq!(chunk.descriptor().chunk(), first_chunk);
    assert_eq!(
        chunk
            .value(crate::LocalRow::new(4))
            .expect("missing residue atom"),
        crate::PropertyValue::Real(0.0, crate::Presence::Unknown)
    );
    assert_eq!(
        chunk
            .value(crate::LocalRow::new(9))
            .expect("non-polymer atom"),
        crate::PropertyValue::Real(0.0, crate::Presence::Inapplicable)
    );
    drop(provider);
    assert!((values[6] - 55.0).abs() < f64::EPSILON);
}
