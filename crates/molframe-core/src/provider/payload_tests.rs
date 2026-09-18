use super::*;

#[test]
fn frame_chunks_share_the_coordinate_pointer_and_own_its_lifetime() {
    let coordinates: CoordinateBlock = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]].into_iter().collect();
    let pointer = coordinates.as_slice().as_ptr();
    let descriptor = ChunkDescriptor::new(
        DatasetId::new(1),
        ChunkId::new(2),
        LogicalRow::new(u64::from(u32::MAX) + 1),
        2,
    )
    .expect("descriptor");
    let chunk = FrameChunk::shared(descriptor, coordinates.clone(), 0..2).expect("frame chunk");

    assert_eq!(chunk.positions().as_ptr(), pointer);
    drop(coordinates);
    assert_eq!(chunk.positions(), &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
}

#[test]
fn local_rows_are_checked_before_global_resolution() {
    let descriptor =
        ChunkDescriptor::new(DatasetId::new(1), ChunkId::new(2), LogicalRow::new(90), 3)
            .expect("descriptor");
    assert_eq!(
        descriptor.resolve(LocalRow::new(2)),
        Ok(LogicalRow::new(92))
    );
    assert!(matches!(
        descriptor.resolve(LocalRow::new(3)),
        Err(ProviderError::LocalRowOutOfRange { .. })
    ));
}

#[test]
fn frame_chunks_reject_descriptor_payload_mismatch() {
    let coordinates: CoordinateBlock = [[1.0, 2.0, 3.0]].into_iter().collect();
    let descriptor =
        ChunkDescriptor::new(DatasetId::new(1), ChunkId::new(2), LogicalRow::new(0), 2)
            .expect("descriptor");
    assert!(matches!(
        FrameChunk::shared(descriptor, coordinates, 0..1),
        Err(ProviderError::PayloadLengthMismatch { .. })
    ));
}
