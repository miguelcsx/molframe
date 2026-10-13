//! Interoperability workloads: graph construction and columnar export.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, structure_with_atoms};
use molframe_interop::{
    AtomTable, EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphOptions, MissingFeaturePolicy,
    NodeFeature, NodeLevel, graph,
};
use molframe_spatial::SpatialBackend;

const GRAPH_TIERS: &[Tier] = &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large];
const TABLE_TIERS: &[Tier] = &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large];

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "interop.graph_radius",
            axis: Axis::Atoms,
            tiers: GRAPH_TIERS,
            run: radius_case,
        },
        Case {
            id: "interop.graph_knn",
            axis: Axis::Atoms,
            tiers: GRAPH_TIERS,
            run: knn_case,
        },
        Case {
            id: "interop.atom_record_batches",
            axis: Axis::Atoms,
            tiers: TABLE_TIERS,
            run: table_case,
        },
    ]
}

fn options(edges: EdgeKind) -> GraphOptions {
    GraphOptions {
        nodes: NodeLevel::Atoms,
        edges,
        direction: EdgeDirection::Undirected,
        node_features: vec![NodeFeature::PositionX, NodeFeature::Element],
        edge_features: vec![EdgeFeature::Distance],
        backend: SpatialBackend::CellList,
        periodic: false,
        missing: MissingFeaturePolicy::Fill(0.0),
    }
}

fn digest(built: &Graph) -> u64 {
    let mut sum = Checksum::default();
    sum.word(built.node_count as u64);
    sum.word(built.edge_count as u64);
    for &index in &built.edge_index {
        sum.word(index.unsigned_abs());
    }
    sum.value()
}

fn graph_case(run: &Run<'_>, atoms: usize, edges: EdgeKind) -> Result<Outcome, String> {
    let structure = structure_with_atoms(atoms as u64);
    let options = options(edges);
    let sample = run.measure(
        || graph(&structure, &options, &run.context).map_err(|e| e.to_string()),
        digest,
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn radius_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    graph_case(run, atoms, EdgeKind::Radius { cutoff: 4.0 })
}

fn knn_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    graph_case(run, atoms, EdgeKind::KNearest { neighbors: 16 })
}

fn table_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let structure = structure_with_atoms(atoms as u64);
    let table = AtomTable::new(&structure);
    let sample = run.measure(
        || table.record_batches().map_err(|e| e.to_string()),
        |batches| {
            let mut sum = Checksum::default();
            for batch in batches {
                sum.word(batch.num_rows() as u64);
                sum.word(batch.get_array_memory_size() as u64);
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}
