//! Criterion coverage for deterministic adapter projections.

use criterion::{Criterion, black_box};
use pdbiox_adapters::TopologyBatch;
use pdbiox_bench::{Sample, structure};
use pdbiox_core::contract::Namespace;
use pdbiox_core::index::ModelIndex;

fn bench_topology_projection(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    c.bench_function("adapters_topology_export/4hhb", |b| {
        b.iter(|| {
            let result = TopologyBatch::from_model(
                black_box(&structure),
                ModelIndex::new(0),
                Namespace::Auth,
            );
            black_box(result)
        });
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_topology_projection(&mut criterion);
    criterion.final_summary();
}
