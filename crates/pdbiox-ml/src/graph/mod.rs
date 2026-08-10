//! Deterministic graph export for tensor graph libraries.

mod build;
mod edges;
mod features;
mod model;
mod nodes;

pub use build::graph;
pub use model::{
    EdgeDirection, EdgeFeature, EdgeKind, FeatureMatrix, Graph, GraphError, GraphOptions,
    MissingFeaturePolicy, NodeFeature, NodeLevel,
};

#[cfg(test)]
#[path = "graph_tests.rs"]
mod tests;
