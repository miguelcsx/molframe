use super::*;
use pdbiox_core::chunk::AtomRecord;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::CoordinateStore;
use pdbiox_core::{
    AltId, ChunkBuilder, Element, OptionalSymbol, Presence, ResidueIndex, StructureData, SymbolId,
};
use pdbiox_query::{Groups, Query};

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.label_entity_id\n_atom_site.pdbx_PDB_model_num\n\
_atom_site.auth_atom_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 1 1 C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1 1 C2 LIG A 1 1 0 0\n\
ATOM 3 C C3 LIG A 1 1 1 C3 LIG A 1 3 0 0\n\
ATOM 4 C C4 LIG A 1 1 1 C4 LIG A 1 0 0 2\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn select(source: &str, backend: SpatialBackend) -> Vec<u32> {
    let context = pdbiox_core::ExecutionContext::default();
    let structure = structure();
    let policy =
        AnalysisPolicy::default().with_identifiers(pdbiox_core::contract::Namespace::Label);
    let resolver = match StructureSpatial::new(&structure, &policy, backend, &context) {
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
        assert_eq!(
            select("within 1.1 of name C1", backend),
            expected,
            "{backend:?}"
        );
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
    let context = pdbiox_core::ExecutionContext::default();
    let structure = structure();
    let policy =
        AnalysisPolicy::default().with_identifiers(pdbiox_core::contract::Namespace::Label);
    let resolver =
        match StructureSpatial::new(&structure, &policy, SpatialBackend::CellList, &context) {
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

#[test]
fn streamed_within_and_iso_match_the_materialized_pair_reference() {
    let context = pdbiox_core::ExecutionContext::default();
    let structure = structure();
    let policy =
        AnalysisPolicy::default().with_identifiers(pdbiox_core::contract::Namespace::Label);
    let universe = AtomSelection::All(structure.atom_count());
    let targets = [
        AtomSelection::Empty,
        AtomSelection::from_sorted(vec![0]),
        AtomSelection::from_sorted(vec![0, 2]),
    ];

    for backend in [SpatialBackend::BruteForce, SpatialBackend::CellList] {
        let resolver = match StructureSpatial::new(&structure, &policy, backend, &context) {
            Ok(resolver) => resolver,
            Err(finding) => panic!("resolver failed for {backend:?}: {finding}"),
        };
        for target in &targets {
            for cutoff in [0.0, 1.1, 3.1] {
                let expected = materialized_within(&resolver, &universe, target, cutoff, backend);
                let actual = match resolver.within(&universe, target, cutoff) {
                    Ok(selection) => selection,
                    Err(finding) => panic!("streamed within failed: {finding}"),
                };
                assert_eq!(actual, expected, "within {backend:?} cutoff {cutoff}");
            }
            for (inner, outer) in [(0.0, 0.0), (0.5, 1.1), (1.1, 3.1)] {
                let expected =
                    materialized_iso(&resolver, &universe, target, inner, outer, backend);
                let actual = match resolver.iso_layer(&universe, target, inner, outer) {
                    Ok(selection) => selection,
                    Err(finding) => panic!("streamed iso failed: {finding}"),
                };
                assert_eq!(
                    actual, expected,
                    "iso {backend:?} interval [{inner}, {outer}]"
                );
            }
        }
    }
}

#[test]
fn sparse_hundred_thousand_atom_reductions_keep_only_linear_state() {
    let context = pdbiox_core::ExecutionContext::default();
    let structure = sparse_paired_structure(50_000);
    let policy = AnalysisPolicy::default();
    let resolver =
        match StructureSpatial::new(&structure, &policy, SpatialBackend::CellList, &context) {
            Ok(resolver) => resolver,
            Err(finding) => panic!("large sparse resolver failed: {finding}"),
        };
    let target = AtomSelection::from_sorted((0_u32..50_000).map(|pair| pair * 2).collect());
    let universe = AtomSelection::from_sorted((0_u32..50_000).map(|pair| pair * 2 + 1).collect());

    let within = match resolver.within(&universe, &target, 0.75) {
        Ok(selection) => selection,
        Err(finding) => panic!("large sparse within failed: {finding}"),
    };
    let iso = match resolver.iso_layer(&universe, &target, 0.4, 0.75) {
        Ok(selection) => selection,
        Err(finding) => panic!("large sparse iso failed: {finding}"),
    };

    assert_eq!(within, universe);
    assert_eq!(iso, universe);
    assert_eq!(resolver.cached_index_count(), 1);
}

fn materialized_within(
    resolver: &StructureSpatial<'_>,
    universe: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
) -> AtomSelection {
    let pairs = match resolver.pairs_with_backend(universe, target, cutoff, backend) {
        Ok(pairs) => pairs,
        Err(finding) => panic!("materialized within reference failed: {finding}"),
    };
    let mut matched = vec![false; resolver.positions().len()];
    for pair in pairs {
        if universe.contains(pair.first) && target.contains(pair.second) {
            matched[test_index(pair.first)] = true;
        }
        if universe.contains(pair.second) && target.contains(pair.first) {
            matched[test_index(pair.second)] = true;
        }
    }
    AtomSelection::from_sorted(
        universe
            .iter()
            .filter(|&atom| matched[test_index(atom)] || target.contains(atom))
            .collect(),
    )
}

fn materialized_iso(
    resolver: &StructureSpatial<'_>,
    universe: &AtomSelection,
    target: &AtomSelection,
    inner: f32,
    outer: f32,
    backend: SpatialBackend,
) -> AtomSelection {
    let pairs = match resolver.pairs_with_backend(universe, target, outer, backend) {
        Ok(pairs) => pairs,
        Err(finding) => panic!("materialized iso reference failed: {finding}"),
    };
    let mut nearest = vec![f32::INFINITY; resolver.positions().len()];
    for pair in pairs {
        if universe.contains(pair.first) && target.contains(pair.second) {
            let index = test_index(pair.first);
            nearest[index] = nearest[index].min(pair.distance_squared);
        }
        if universe.contains(pair.second) && target.contains(pair.first) {
            let index = test_index(pair.second);
            nearest[index] = nearest[index].min(pair.distance_squared);
        }
    }
    let inner_squared = inner * inner;
    let outer_squared = outer * outer;
    AtomSelection::from_sorted(
        universe
            .iter()
            .filter(|&atom| {
                let squared = if target.contains(atom) {
                    0.0
                } else {
                    nearest[test_index(atom)]
                };
                squared >= inner_squared && squared <= outer_squared
            })
            .collect(),
    )
}

fn sparse_paired_structure(pair_count: u32) -> Structure {
    let atom_count = pair_count
        .checked_mul(2)
        .unwrap_or_else(|| panic!("test atom count overflow"));
    let mut builder = ChunkBuilder::new();
    builder.reserve(test_index(atom_count));
    builder.start_model(0);
    let mut base = 0.0_f32;
    for pair in 0..pair_count {
        let first = pair * 2;
        builder.push(test_atom([base, 0.0, 0.0], first));
        builder.push(test_atom([base + 0.5, 0.0, 0.0], first + 1));
        base += 4.0;
    }
    let (chunks, coords) = builder.finish();
    let mut data = StructureData::empty();
    data.chunks = chunks.into();
    data.coords = CoordinateStore::Single(coords);
    Structure::new(data)
}

fn test_atom(position: [f32; 3], atom: u32) -> AtomRecord {
    AtomRecord {
        position: Some(position),
        element: Element::CARBON,
        atom_name: SymbolId::from_raw(0),
        auth_atom_name: OptionalSymbol::NONE,
        alternate_component_id: OptionalSymbol::NONE,
        alt_id: AltId::BLANK,
        residue: ResidueIndex::new(0),
        occupancy: (1.0, Presence::Present),
        b_factor: (0.0, Presence::Present),
        formal_charge: (0, Presence::Inapplicable),
        atom_site_id: atom + 1,
    }
}

fn test_index(atom: u32) -> usize {
    match usize::try_from(atom) {
        Ok(index) => index,
        Err(error) => panic!("test atom index exceeds usize: {error}"),
    }
}
