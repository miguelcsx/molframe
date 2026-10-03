use crate::view::tests::{ENTRY, attached};
use crate::{AssemblyExt, AssemblyNeighbor};
use molframe_core::{ExecutionContext, ModelIndex};
use molframe_spatial::SpatialBackend;

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
    match view.neighbors(
        ModelIndex::new(0),
        3.0,
        backend,
        &ExecutionContext::default(),
    ) {
        Ok(neighbors) => neighbors,
        Err(finding) => panic!("spatial search failed: {finding}"),
    }
}

fn boundary_entry(carbon_x: f32) -> String {
    ENTRY
        .replace(
            "1 C CA GLY A 1 1 X 1 0 0",
            &format!("1 C C GLY A 1 1 X {carbon_x} 0 0"),
        )
        .replace("2 C CA GLY B 1 1 Y 0 2 0", "2 N N GLY B 1 1 Y 0 0 0")
        .replace("R 0 -1 0 0 1 0 0 0 0 0 1 0", "R 1 0 0 0 0 1 0 0 0 0 1 0")
        .replace("1 '(T)(R)' 'A,B'", "1 R A\n1 T B")
}

#[test]
fn a_covalent_pair_crossing_a_ten_angstrom_placement_boundary_is_found() {
    let structure = attached(&boundary_entry(8.67));
    let view = structure.assembly("1").expect("assembly exists");
    let links = view
        .covalent_links(ModelIndex::new(0), &ExecutionContext::default())
        .expect("search succeeds");
    assert_eq!(
        links,
        [crate::AssemblyBond {
            first_instance: molframe_core::InstanceId::new(0),
            first_atom: molframe_core::AtomIndex::new(0),
            second_instance: molframe_core::InstanceId::new(1),
            second_atom: molframe_core::AtomIndex::new(1),
            order: molframe_core::BondOrder::Single,
        }]
    );
}

#[test]
fn a_three_point_two_angstrom_contact_is_not_a_covalent_link() {
    let structure = attached(&boundary_entry(6.8));
    let view = structure.assembly("1").expect("assembly exists");
    assert!(
        view.covalent_links(ModelIndex::new(0), &ExecutionContext::default())
            .expect("search succeeds")
            .is_empty()
    );
}

#[test]
fn same_instance_pairs_are_not_cross_instance_links() {
    let entry = boundary_entry(8.67)
        .replace("2 N N GLY B 1 1 Y 0 0 0", "2 N N GLY A 2 2 X 10 0 0")
        .replace("1 T B", "");
    let structure = attached(&entry);
    let view = structure.assembly("1").expect("assembly exists");
    assert!(
        view.covalent_links(ModelIndex::new(0), &ExecutionContext::default())
            .expect("search succeeds")
            .is_empty()
    );
}

#[test]
fn deposited_4hhb_assembly_has_no_cross_instance_covalent_links() {
    let bytes = molframe_bench::Sample::Medium
        .cif()
        .expect("committed CIF fixture");
    let structure = attached(std::str::from_utf8(bytes).expect("fixture is UTF-8"));
    let view = structure.assembly("1").expect("deposited assembly exists");
    assert!(
        view.covalent_links(ModelIndex::new(0), &ExecutionContext::default())
            .expect("valid real assembly search")
            .is_empty()
    );
}
