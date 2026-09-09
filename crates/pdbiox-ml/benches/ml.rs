//! Criterion coverage for Arrow table and graph materialisation.

use arrow::ffi_stream::ArrowArrayStreamReader;
use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, structure};
use pdbiox_core::execution::ExecutionContext;
use pdbiox_ml::{
    AtomTable, EdgeDirection, EdgeFeature, EdgeKind, GraphOptions, MissingFeaturePolicy,
    NodeFeature, NodeLevel, graph,
};
use pdbiox_spatial::SpatialBackend;

fn bench_arrow(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let table = AtomTable::new(&structure);
    c.bench_function("ml_atom_record_batches/4hhb", |b| {
        b.iter(|| black_box(table.record_batches()));
    });
    c.bench_function("ml_atom_arrow_stream/4hhb", |b| {
        b.iter(|| black_box(table.arrow_stream()));
    });
    c.bench_function("ml_atom_arrow_stream_consume/4hhb", |b| {
        b.iter(|| {
            let stream = match table.arrow_stream() {
                Ok(stream) => stream,
                Err(error) => panic!("Arrow stream creation failed: {error}"),
            };
            let reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
                Ok(reader) => reader,
                Err(error) => panic!("Arrow stream import failed: {error}"),
            };
            let mut rows = 0_usize;
            for batch in reader {
                rows += match batch {
                    Ok(batch) => batch.num_rows(),
                    Err(error) => panic!("Arrow stream pull failed: {error}"),
                };
            }
            black_box(rows)
        });
    });
}

fn bench_graph(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let structure = structure(Sample::Small);
    let options = GraphOptions {
        nodes: NodeLevel::Atoms,
        edges: EdgeKind::Radius { cutoff: 4.0 },
        direction: EdgeDirection::Undirected,
        node_features: vec![NodeFeature::PositionX, NodeFeature::Element],
        edge_features: vec![EdgeFeature::Distance],
        backend: SpatialBackend::CellList,
        periodic: false,
        missing: MissingFeaturePolicy::Fill(0.0),
    };
    if let Err(error) = graph(&structure, &options, &context) {
        panic!("ML graph benchmark setup failed: {error}");
    }
    c.bench_function("ml_graph/1ubq", |b| {
        b.iter(|| black_box(graph(&structure, &options, &context)));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_arrow(&mut criterion);
    bench_graph(&mut criterion);
    criterion.final_summary();
}
