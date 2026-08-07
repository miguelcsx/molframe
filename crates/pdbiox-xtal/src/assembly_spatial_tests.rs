use crate::view::tests::{ENTRY, attached};
use crate::{AssemblyExt, AssemblyNeighbor};
use pdbiox_core::ModelIndex;
use pdbiox_spatial::SpatialBackend;

#[test]
fn transformed_assembly_pairs_are_identical_across_every_backend() {
    let structure = attached(ENTRY);
    let view = match structure.assembly("1") {
        Ok(view) => view,
        Err(finding) => panic!("assembly failed: {finding}"),
    };
    let expected = neighbors(&view, SpatialBackend::BruteForce);
    assert_eq!(expected.len(), 1);
    assert!((expected[0].distance_squared - 5.0).abs() < 1e-6);
    for backend in [
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
        SpatialBackend::Auto,
    ] {
        assert_eq!(neighbors(&view, backend), expected);
    }
}

fn neighbors(view: &crate::AssemblyView, backend: SpatialBackend) -> Vec<AssemblyNeighbor> {
    match view.neighbors(ModelIndex::new(0), 3.0, backend) {
        Ok(neighbors) => neighbors,
        Err(finding) => panic!("spatial search failed: {finding}"),
    }
}
