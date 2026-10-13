//! Criterion coverage for selection over a million-atom structure.

use criterion::{Criterion, black_box};
use molframe_bench::structure_with_atoms;
use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_query::{Groups, Query};
use molframe_spatial::{SpatialBackend, StructureSpatial};

const ATOMS: u64 = 1_000_000;

const QUERIES: [(&str, &str); 8] = [
    ("name_ca", "name CA"),
    ("and_chain", "name CA and chain A"),
    ("numeric", "bfactor > 20 and occupancy > 0.5"),
    ("element_zn", "element ZN"),
    ("not_or", "not (name CA or name N) or resname HOH"),
    ("within", "within 5 of resname ZN"),
    ("byres", "byres (name CA and bfactor > 40)"),
    ("byres_sparse", "byres element ZN"),
];

fn bench_scale(c: &mut Criterion) {
    let structure = structure_with_atoms(ATOMS);
    let policy = AnalysisPolicy::default();
    let groups = Groups::new();
    let context = ExecutionContext::default();
    let spatial =
        match StructureSpatial::new(&structure, &policy, SpatialBackend::CellList, &context) {
            Ok(spatial) => spatial,
            Err(error) => panic!("scale spatial fixture failed: {error:?}"),
        };
    let mut group = c.benchmark_group("query_scale_1m");
    group.sample_size(10);
    for (label, text) in QUERIES {
        let query = match Query::compile(text) {
            Ok(query) => query,
            Err(findings) => panic!("scale query {text} failed: {findings:?}"),
        };
        group.bench_function(label, |b| {
            b.iter(|| black_box(query.evaluate(&structure, &policy, &groups, Some(&spatial))));
        });
    }
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_scale(&mut criterion);
    criterion.final_summary();
}
