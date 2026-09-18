use super::*;
use crate::provider::LocalRow;

#[test]
fn a_terabyte_scale_catalog_does_not_materialise_chunk_metadata() {
    let logical_rows = 1_u64 << 40;
    let dataset = DatasetDescriptor::regular(
        DatasetId::new(9),
        PayloadKind::Frame,
        logical_rows,
        1,
        ChunkId::new(1_u64 << 48),
    )
    .expect("regular layout is representable");
    let catalog = DatasetCatalog::new(vec![dataset]).expect("one dataset is valid");

    assert_eq!(catalog.len(), 1);
    assert!(dataset.chunk_count() > u64::from(u32::MAX));
    assert_eq!(catalog.datasets().len(), 1);
    assert_eq!(
        std::mem::size_of_val(&catalog),
        std::mem::size_of::<DatasetCatalog>()
    );
}

#[test]
fn regular_chunks_resolve_rows_beyond_u32_without_saturation() {
    let dataset = DatasetDescriptor::regular(
        DatasetId::new(u64::MAX - 3),
        PayloadKind::Property,
        u64::from(u32::MAX) + 10_000,
        4096,
        ChunkId::new(70_000),
    )
    .expect("layout is representable");
    let ordinal = u64::from(u32::MAX) / 4096 + 1;
    let chunk = ChunkId::new(dataset.first_chunk().get() + ordinal);
    let descriptor = dataset.regular_chunk(chunk).expect("chunk exists");

    assert!(descriptor.logical_start().get() > u64::from(u32::MAX) - 4096);
    let global = descriptor
        .resolve(LocalRow::new(descriptor.rows() - 1))
        .expect("tail local row resolves");
    assert!(global.get() > u64::from(u32::MAX));
}

#[test]
fn identity_overflow_is_an_error_instead_of_a_saturated_range() {
    let result = DatasetDescriptor::regular(
        DatasetId::new(1),
        PayloadKind::Structure,
        2,
        1,
        ChunkId::new(u64::MAX),
    );
    assert!(matches!(
        result,
        Err(ProviderError::IdentityOverflow {
            field: "chunk identity",
            ..
        })
    ));
}

#[test]
fn catalogs_reject_overlapping_chunk_namespaces() {
    let first = DatasetDescriptor::regular(
        DatasetId::new(1),
        PayloadKind::Structure,
        10,
        5,
        ChunkId::new(10),
    )
    .expect("first layout");
    let second = DatasetDescriptor::regular(
        DatasetId::new(2),
        PayloadKind::Frame,
        10,
        5,
        ChunkId::new(11),
    )
    .expect("second layout");
    assert!(matches!(
        DatasetCatalog::new(vec![first, second]),
        Err(ProviderError::OverlappingChunkRange { .. })
    ));
}
