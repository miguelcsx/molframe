use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_query::{Groups, Query};

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1 0 0\n\
ATOM 3 C C3 LIG A 1 3 0 0\n\
ATOM 4 C C4 LIG A 1 0 0 2\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn select(source: &str, backend: SpatialBackend) -> Vec<u32> {
    let structure = structure();
    let policy = AnalysisPolicy::default();
    let resolver = match StructureSpatial::new(&structure, &policy, backend) {
        Ok(resolver) => resolver,
        Err(finding) => panic!("resolver failed: {finding}"),
    };
    let query = match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    let result = query.evaluate(&structure, &policy, &Groups::new(), Some(&resolver));
    match result {
        Ok(result) => result.selection.iter().collect(),
        Err(findings) => panic!("evaluation failed: {findings:?}"),
    }
}

#[test]
fn textual_geometric_queries_are_identical_across_backends() {
    let expected = select("within 1.1 of name C1", SpatialBackend::BruteForce);
    assert_eq!(expected, vec![0, 1]);
    for backend in [
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        assert_eq!(select("within 1.1 of name C1", backend), expected);
    }
}

#[test]
fn around_beyond_point_spheres_shells_and_cylinders_have_distinct_semantics() {
    assert_eq!(
        select("around 1.1 name C1", SpatialBackend::CellList),
        vec![1]
    );
    assert_eq!(
        select("beyond 1.1 of name C1", SpatialBackend::CellList),
        vec![2, 3]
    );
    assert_eq!(
        select("point 0 0 0 1.1", SpatialBackend::CellList),
        vec![0, 1]
    );
    assert_eq!(
        select("sphzone 1.1 name C1", SpatialBackend::CellList),
        vec![0, 1]
    );
    assert_eq!(
        select("sphlayer 0.5 1.1 name C1", SpatialBackend::CellList),
        vec![1]
    );
    assert_eq!(
        select("isolayer 0.5 1.1 name C1", SpatialBackend::CellList),
        vec![1]
    );
    assert_eq!(
        select("cyzone 1.1 0.5 -0.5 name C1", SpatialBackend::CellList),
        vec![0, 1]
    );
    assert_eq!(
        select("cylayer 0.5 1.1 0.5 -0.5 name C1", SpatialBackend::CellList),
        vec![1]
    );
}

#[test]
fn repeated_structure_queries_reuse_the_same_generation_bound_index() {
    let structure = structure();
    let policy = AnalysisPolicy::default();
    let resolver = match StructureSpatial::new(&structure, &policy, SpatialBackend::CellList) {
        Ok(resolver) => resolver,
        Err(finding) => panic!("resolver failed: {finding}"),
    };
    let query = match Query::compile("within 1.1 of name C1") {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    for _ in 0..2 {
        if let Err(findings) = query.evaluate(&structure, &policy, &Groups::new(), Some(&resolver))
        {
            panic!("evaluation failed: {findings:?}")
        }
    }
    assert_eq!(resolver.cached_index_count(), 1);
}
