use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_spatial::SpatialBackend;

const SOURCE: &str = "data_graph\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 N N GLY A 1 1 0 0\n\
ATOM 3 O O HOH B 2 4 0 0\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n1 covale A 1 GLY CA A 1 GLY N SING\n";

fn structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn options(edges: EdgeKind) -> GraphOptions {
    GraphOptions {
        nodes: NodeLevel::Atoms,
        edges,
        direction: EdgeDirection::Symmetric,
        node_features: vec![
            NodeFeature::Element,
            NodeFeature::PositionX,
            NodeFeature::PositionY,
            NodeFeature::PositionZ,
        ],
        edge_features: vec![EdgeFeature::Distance],
        missing: MissingFeaturePolicy::Error,
        backend: SpatialBackend::Auto,
        periodic: false,
    }
}

#[test]
fn bond_graph_matches_pyg_shapes_and_typed_features() {
    let mut options = options(EdgeKind::Bonds);
    options.edge_features.push(EdgeFeature::BondOrder);
    let graph = graph(&structure(), &options).expect("bond graph");
    assert_eq!((graph.node_count, graph.edge_count), (3, 2));
    assert_eq!(&*graph.edge_index, &[0, 1, 1, 0]);
    assert_eq!(
        (graph.node_features.rows, graph.node_features.columns),
        (3, 4)
    );
    assert_eq!(
        (graph.edge_features.rows, graph.edge_features.columns),
        (2, 2)
    );
    assert_eq!(&*graph.edge_features.values, &[1.0, 1.0, 1.0, 1.0]);
    assert_eq!(graph.cost(), crate::ExportCost::Copy);
}

#[test]
fn radius_and_nearest_graphs_are_deterministic() {
    let radius = options(EdgeKind::Radius { cutoff: 4.1 });
    let first = graph(&structure(), &radius).expect("first radius graph");
    let second = graph(&structure(), &radius).expect("second radius graph");
    assert_eq!(first, second);

    let nearest = options(EdgeKind::KNearest { neighbors: 1 });
    let first = graph(&structure(), &nearest).expect("first nearest graph");
    let second = graph(&structure(), &nearest).expect("second nearest graph");
    assert_eq!(first, second);
}

#[test]
fn missing_features_require_an_explicit_finite_fill() {
    let mut strict = options(EdgeKind::Radius { cutoff: 2.0 });
    strict.node_features = vec![NodeFeature::PartialCharge];
    assert!(matches!(
        graph(&structure(), &strict),
        Err(GraphError::MissingFeature { .. })
    ));
    strict.missing = MissingFeaturePolicy::Fill(-7.0);
    let filled = graph(&structure(), &strict).expect("explicit fill graph");
    assert_eq!(&*filled.node_features.values, &[-7.0; 3]);
}

#[test]
fn residue_graph_aggregates_positions_without_atom_level_features() {
    let mut options = options(EdgeKind::Radius { cutoff: 5.0 });
    options.nodes = NodeLevel::Residues;
    options.node_features = vec![NodeFeature::AtomCount, NodeFeature::PositionX];
    let graph = graph(&structure(), &options).expect("residue graph");
    assert_eq!(graph.node_count, 2);
    assert_eq!(&*graph.node_features.values, &[2.0, 0.5, 1.0, 4.0]);
}
