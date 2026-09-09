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
