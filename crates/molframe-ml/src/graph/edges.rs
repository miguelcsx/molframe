//! Spatial and topological edge construction with stable ordering.

use super::nodes::Nodes;
use super::{EdgeDirection, EdgeKind, GraphError, GraphOptions, NodeLevel};
use molframe_core::{AtomSelection, BondOrder, ExecutionContext, Structure};
use molframe_spatial::{KdTree, NeighborPair, PeriodicBox, SpatialBackend, pairs_within};
use std::collections::{BTreeMap, BTreeSet};

const BONDS_FEATURE: &str = "bonds";
const DIRECTED_GEOMETRY_FEATURE: &str = "directed_non_knn_edges";
const RESIDUE_BOND_FEATURE: &str = "bond_edges_for_residues";

#[derive(Clone, Copy)]
pub(super) struct Edge {
    pub(super) source: u32,
    pub(super) target: u32,
    pub(super) distance_squared: f32,
    pub(super) bond_order: Option<BondOrder>,
}

pub(super) fn build(
    structure: &Structure,
    nodes: &Nodes,
    options: &GraphOptions,
    context: &ExecutionContext,
) -> Result<Vec<Edge>, GraphError> {
    let periodic = periodic_box(structure, options.periodic)?;
    let base = match options.edges {
        EdgeKind::Bonds => bonds(structure, nodes, options.nodes, periodic.as_ref())?,
        EdgeKind::Contacts { cutoff } => geometric(
            structure,
            nodes,
            cutoff,
            options.backend,
            periodic.as_ref(),
            true,
            context,
        )?,
        EdgeKind::Radius { cutoff } => geometric(
            structure,
            nodes,
            cutoff,
            options.backend,
            periodic.as_ref(),
            false,
            context,
        )?,
        EdgeKind::KNearest { neighbors } => nearest(nodes, neighbors, periodic.as_ref())?,
    };
    orient(base, options.direction, options.edges)
}

fn bonds(
    structure: &Structure,
    nodes: &Nodes,
    level: NodeLevel,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<Edge>, GraphError> {
    if level != NodeLevel::Atoms {
        return Err(GraphError::UnsupportedFeature {
            feature: RESIDUE_BOND_FEATURE,
        });
    }
    if !structure.data().bonds.is_available() {
        return Err(GraphError::MissingFeature {
            feature: BONDS_FEATURE,
            row: 0,
        });
    }
    structure
        .data()
        .bonds
        .iter()
        .map(|bond| {
            edge(
                nodes,
                bond.atom_a.get(),
                bond.atom_b.get(),
                periodic,
                Some(bond.order),
            )
        })
        .collect()
}

fn geometric(
    structure: &Structure,
    nodes: &Nodes,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    exclude_bonds: bool,
    context: &ExecutionContext,
) -> Result<Vec<Edge>, GraphError> {
    if exclude_bonds && !structure.data().bonds.is_available() {
        return Err(GraphError::MissingFeature {
            feature: BONDS_FEATURE,
            row: 0,
        });
    }
    let count = u32::try_from(nodes.positions.len()).map_err(|_| GraphError::IndexOverflow)?;
    let selection = AtomSelection::All(count);
    let pairs = pairs_within(
        &nodes.positions,
        &selection,
        &selection,
        cutoff,
        backend,
        periodic,
        context,
    )?;
    let excluded = exclude_bonds.then(|| bonded_nodes(structure, nodes));
    Ok(pairs
        .into_iter()
        .filter(|pair| {
            excluded
                .as_ref()
                .is_none_or(|bonds| !bonds.contains(&(pair.first, pair.second)))
        })
        .map(pair_edge)
        .collect())
}

fn nearest(
    nodes: &Nodes,
    neighbors: usize,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<Edge>, GraphError> {
    if neighbors == 0 {
        return Err(GraphError::InvalidParameter);
    }
    let targets: Vec<_> = (0..nodes.positions.len())
        .map(|node| u32::try_from(node).map_err(|_| GraphError::IndexOverflow))
        .collect::<Result<_, _>>()?;
    let tree = KdTree::build(&nodes.positions, &targets, periodic.copied())?;
    let capacity = nodes
        .positions
        .len()
        .checked_mul(neighbors)
        .ok_or(GraphError::ResourceLimit)?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(capacity)
        .map_err(|_| GraphError::ResourceLimit)?;
    for source in targets {
        for (target, distance_squared) in tree.k_nearest(source, neighbors)? {
            edges.push(Edge {
                source,
                target,
                distance_squared,
                bond_order: None,
            });
        }
    }
    edges.sort_unstable_by_key(|edge| (edge.source, edge.target));
    edges.dedup_by_key(|edge| (edge.source, edge.target));
    Ok(edges)
}

fn orient(
    edges: Vec<Edge>,
    direction: EdgeDirection,
    kind: EdgeKind,
) -> Result<Vec<Edge>, GraphError> {
    if direction == EdgeDirection::Directed && !matches!(kind, EdgeKind::KNearest { .. }) {
        return Err(GraphError::UnsupportedFeature {
            feature: DIRECTED_GEOMETRY_FEATURE,
        });
    }
    let mut output = BTreeMap::new();
    for edge in edges {
        let canonical = canonical_edge(edge);
        match direction {
            EdgeDirection::Undirected => {
                output.insert((canonical.source, canonical.target), canonical);
            }
            EdgeDirection::Symmetric => {
                output.insert((canonical.source, canonical.target), canonical);
                let reverse = Edge {
                    source: canonical.target,
                    target: canonical.source,
                    ..canonical
                };
                output.insert((reverse.source, reverse.target), reverse);
            }
            EdgeDirection::Directed => {
                output.insert((edge.source, edge.target), edge);
            }
        }
    }
    Ok(output.into_values().collect())
}

fn edge(
    nodes: &Nodes,
    source: u32,
    target: u32,
    periodic: Option<&PeriodicBox>,
    bond_order: Option<BondOrder>,
) -> Result<Edge, GraphError> {
    let Some(left) = nodes.positions.get(source as usize).copied() else {
        return Err(GraphError::IndexOverflow);
    };
    let Some(right) = nodes.positions.get(target as usize).copied() else {
        return Err(GraphError::IndexOverflow);
    };
    let distance_squared = match periodic {
        Some(periodic) => periodic.distance_squared(left, right),
        None => left
            .iter()
            .zip(right)
            .map(|(left, right)| (left - right).powi(2))
            .sum(),
    };
    Ok(Edge {
        source,
        target,
        distance_squared,
        bond_order,
    })
}

fn pair_edge(pair: NeighborPair) -> Edge {
    Edge {
        source: pair.first,
        target: pair.second,
        distance_squared: pair.distance_squared,
        bond_order: None,
    }
}

fn canonical_edge(edge: Edge) -> Edge {
    if edge.source <= edge.target {
        edge
    } else {
        Edge {
            source: edge.target,
            target: edge.source,
            ..edge
        }
    }
}

fn bonded_nodes(structure: &Structure, nodes: &Nodes) -> BTreeSet<(u32, u32)> {
    structure
        .data()
        .bonds
        .iter()
        .filter_map(|bond| {
            let left = *nodes.atom_to_node.get(bond.atom_a.as_usize())?;
            let right = *nodes.atom_to_node.get(bond.atom_b.as_usize())?;
            (left != right).then_some(if left < right {
                (left, right)
            } else {
                (right, left)
            })
        })
        .collect()
}

fn periodic_box(structure: &Structure, periodic: bool) -> Result<Option<PeriodicBox>, GraphError> {
    if !periodic {
        return Ok(None);
    }
    let cell = structure.data().cell.ok_or(GraphError::MissingCell)?;
    PeriodicBox::from_cell(cell)
        .map(Some)
        .map_err(|_| GraphError::MissingCell)
}
