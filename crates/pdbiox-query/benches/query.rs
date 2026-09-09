//! Criterion coverage for reusable typed and textual selections.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, structure};
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_query::{Groups, Query, col};
use pdbiox_spatial::{SpatialBackend, StructureSpatial};

fn bench_queries(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let structure = structure(Sample::Medium);
    let policy = AnalysisPolicy::default();
    let groups = Groups::new();
    let spatial =
        match StructureSpatial::new(&structure, &policy, SpatialBackend::CellList, &context) {
            Ok(spatial) => spatial,
            Err(error) => panic!("query spatial benchmark fixture failed: {error:?}"),
        };
    let typed = Query::from_builder(col::bfactor().gt(20.0));
    let geometric = Query::from_builder(col::within(5.0, col::chain().eq("A")));
    let residue_expansion = Query::from_builder(col::by_residue(col::name().eq("CA")));
    let textual = match Query::compile("chain A and bfactor > 20") {
        Ok(query) => query,
        Err(findings) => panic!("query bench fixture failed: {findings:?}"),
    };
    for query in [&typed, &textual, &geometric, &residue_expansion] {
        if let Err(findings) = query.evaluate(&structure, &policy, &groups, Some(&spatial)) {
            panic!("query benchmark evaluation failed: {findings:?}");
        }
    }
    let mut group = c.benchmark_group("query_evaluation");
    group.bench_function("typed", |b| {
        b.iter(|| black_box(typed.evaluate(&structure, &policy, &groups, Some(&spatial))));
    });
    group.bench_function("textual", |b| {
        b.iter(|| black_box(textual.evaluate(&structure, &policy, &groups, Some(&spatial))));
    });
    group.bench_function("geometric_within", |b| {
        b.iter(|| black_box(geometric.evaluate(&structure, &policy, &groups, Some(&spatial))));
    });
    group.bench_function("residue_expansion", |b| {
        b.iter(|| {
            black_box(residue_expansion.evaluate(&structure, &policy, &groups, Some(&spatial)))
        });
    });
    group.finish();

    c.bench_function("query_compile/complex", |b| {
        b.iter(|| {
            black_box(Query::compile(
                "protein and name C* and bfactor > 20 and within 5 of (resname HEM)",
            ))
        });
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_queries(&mut criterion);
    criterion.final_summary();
}
