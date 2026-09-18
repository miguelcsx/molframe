//! Linear assembly of framework-neutral graph tensors.

use super::features;
use super::nodes;
use super::{EdgeFeature, EdgeKind, Graph, GraphError, GraphOptions, NodeFeature, NodeLevel};
use molframe_core::{ExecutionContext, Structure};
use std::collections::BTreeSet;

/// Materialises graph indices and selected feature matrices.
///
/// Fixed-radius construction delegates to the spatial planner and runs in
/// `O(nodes + edges)` for its cell-list path. Feature materialisation is one
/// row-major pass over the requested columns. Output is a copy because graph
/// tensor layout differs from structure topology and encoded columns.
///
/// # Errors
///
/// Returns an error for contradictory options, absent requested features,
/// unusable periodic cells, invalid geometric parameters or index overflow.
pub fn graph(
    structure: &Structure,
    options: &GraphOptions,
    context: &ExecutionContext,
) -> Result<Graph, GraphError> {
    validate(options)?;
    let nodes = nodes::project(structure, options.nodes)?;
    let edges = super::edges::build(structure, &nodes, options, context)?;
    let node_features = features::nodes(
        structure,
        &nodes,
        options.nodes,
        &options.node_features,
        options.missing,
    )?;
    let edge_features = features::edges(&edges, &options.edge_features, options.missing)?;
    let edge_count = edges.len();
    let capacity = edge_count.checked_mul(2).ok_or(GraphError::ResourceLimit)?;
    let mut edge_index = Vec::new();
    edge_index
        .try_reserve_exact(capacity)
        .map_err(|_| GraphError::ResourceLimit)?;
    edge_index.extend(edges.iter().map(|edge| i64::from(edge.source)));
    edge_index.extend(edges.iter().map(|edge| i64::from(edge.target)));
    Ok(Graph {
        node_count: nodes.positions.len(),
        edge_index: edge_index.into_boxed_slice(),
        edge_count,
        node_features,
        edge_features,
    })
}

fn validate(options: &GraphOptions) -> Result<(), GraphError> {
    unique(&options.node_features, features::node_label)?;
    unique(&options.edge_features, features::edge_label)?;
    match options.edges {
        EdgeKind::Contacts { cutoff } | EdgeKind::Radius { cutoff }
            if !cutoff.is_finite() || cutoff <= 0.0 =>
        {
            return Err(GraphError::InvalidParameter);
        }
        EdgeKind::KNearest { neighbors: 0 } => return Err(GraphError::InvalidParameter),
        _ => {}
    }
    if options.nodes == NodeLevel::Residues && options.node_features.contains(&NodeFeature::Element)
    {
        return Err(GraphError::UnsupportedFeature { feature: "element" });
    }
    if options.nodes == NodeLevel::Atoms && options.node_features.contains(&NodeFeature::AtomCount)
    {
        return Err(GraphError::UnsupportedFeature {
            feature: "atom_count",
        });
    }
    if options.edge_features.contains(&EdgeFeature::BondOrder)
        && !matches!(options.edges, EdgeKind::Bonds)
    {
        return Err(GraphError::UnsupportedFeature {
            feature: "bond_order",
        });
    }
    Ok(())
}

fn unique<T: Copy + Ord>(
    values: &[T],
    label: impl Fn(T) -> &'static str,
) -> Result<(), GraphError> {
    let mut seen = BTreeSet::new();
    for value in values {
        if !seen.insert(*value) {
            return Err(GraphError::UnsupportedFeature {
                feature: label(*value),
            });
        }
    }
    Ok(())
}
