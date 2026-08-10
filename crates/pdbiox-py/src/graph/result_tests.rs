use super::*;
use crate::graph::options::PyGraphOptions;
use crate::structure::PyStructure;
use numpy::PyUntypedArrayMethods;
use pdbiox::{
    EdgeDirection, EdgeFeature, EdgeKind, GraphOptions, MissingFeaturePolicy, NodeFeature,
    NodeLevel, SpatialBackend,
};
use std::path::PathBuf;

#[test]
fn graph_binding_materialises_native_arrays_once_with_typed_schema() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match pdbiox::read(path) {
        Ok(structure) => structure,
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    };
    let options = PyGraphOptions {
        inner: GraphOptions {
            nodes: NodeLevel::Atoms,
            edges: EdgeKind::Radius { cutoff: 4.0 },
            direction: EdgeDirection::Symmetric,
            node_features: vec![NodeFeature::Element],
            edge_features: vec![EdgeFeature::Distance],
            missing: MissingFeaturePolicy::Error,
            backend: SpatialBackend::BruteForce,
            periodic: false,
        },
    };
    let structure = PyStructure::new(structure);
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        let graph = match structure.graph(py, &options) {
            Ok(graph) => graph,
            Err(error) => panic!("native graph export failed: {error}"),
        };
        assert_eq!(graph.edge_index.bind(py).shape(), [2, graph.edge_count]);
        assert_eq!(graph.node_features.bind(py).shape(), [graph.node_count, 1]);
        assert_eq!(graph.edge_features.bind(py).shape(), [graph.edge_count, 1]);
        assert_eq!(graph.node_feature_schema, vec![PyNodeFeature::Element]);
        assert_eq!(graph.edge_feature_schema, vec![PyEdgeFeature::Distance]);
    });
}
