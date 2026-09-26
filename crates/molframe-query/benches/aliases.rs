//! Named-query resolution, fingerprinting, and compile-once reuse at scale.
//!
//! The "uncached" rows compile the query text on every evaluation, which is
//! what a consumer that only holds query text pays without a producer-side
//! cache; the "cached" rows evaluate a query compiled once. The sparse rows run
//! over 4HHB and over a synthetic structure of about a million atoms tiled from
//! 1AON, selecting one atom name of one residue number of one chain.

use criterion::{BenchmarkId, Criterion, black_box};
use molframe_bench::{Sample, Seed, SyntheticCifSource, Tile, input, structure};
use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::ReadOptions;
use molframe_core::structure::Structure;
use molframe_query::{Groups, Query, QueryAliases};
use molframe_spatial::{SpatialBackend, StructureSpatial};
use std::io::Read as _;

fn compile(source: &str) -> Query {
    match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("bench query {source:?} failed: {findings:?}"),
    }
}

fn aliases(depth: usize) -> (QueryAliases, Query) {
    let mut aliases = QueryAliases::new();
    let _ = aliases.define("n0", compile("resname HEM"));
    for level in 1..depth {
        let _ = aliases.define(
            format!("n{level}"),
            compile(&format!("byres (within 4 of $n{}) and protein", level - 1)),
        );
    }
    (aliases, compile(&format!("$n{}", depth - 1)))
}

fn synthetic(atoms: u64) -> Structure {
    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let mut source = SyntheticCifSource::with_atoms(tile, atoms);
    let mut bytes = Vec::new();
    if let Err(error) = source.read_to_end(&mut bytes) {
        panic!("synthetic stream failed: {error}");
    }
    match molframe_cif::read(&input(&bytes), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("synthetic structure failed: {findings:?}"),
    }
}

fn bench_language(c: &mut Criterion) {
    let source = "byres (within 5 of (resname HEM and chain A)) and protein and not water";
    let query = compile(source);
    c.bench_function("query_language/compile", |b| {
        b.iter(|| black_box(Query::compile(black_box(source))));
    });
    c.bench_function("query_language/fingerprint", |b| {
        b.iter(|| black_box(query.fingerprint()));
    });
    let mut group = c.benchmark_group("query_aliases/resolve");
    for depth in [1_usize, 4, 16] {
        let (aliases, query) = aliases(depth);
        group.bench_with_input(BenchmarkId::from_parameter(depth), &depth, |b, _| {
            b.iter(|| black_box(aliases.resolve(black_box(&query))));
        });
    }
    group.finish();
}

fn bench_reuse(c: &mut Criterion) {
    let policy = AnalysisPolicy::default();
    let groups = Groups::new();
    let context = ExecutionContext::default();
    let mut group = c.benchmark_group("query_reuse");
    group.sample_size(20);
    for (label, structure) in [
        ("4hhb", structure(Sample::Medium)),
        ("synthetic_1m", synthetic(1_000_000)),
    ] {
        let spatial =
            match StructureSpatial::new(&structure, &policy, SpatialBackend::Auto, &context) {
                Ok(spatial) => spatial,
                Err(error) => panic!("spatial fixture failed: {error:?}"),
            };
        let source = "chain A and resid 10 and name CA";
        let compiled = compile(source);
        group.bench_function(BenchmarkId::new("uncached_sparse", label), |b| {
            b.iter(|| {
                let query = compile(black_box(source));
                black_box(query.evaluate(&structure, &policy, &groups, Some(&spatial)))
            });
        });
        group.bench_function(BenchmarkId::new("cached_sparse", label), |b| {
            b.iter(|| black_box(compiled.evaluate(&structure, &policy, &groups, Some(&spatial))));
        });
        let physical = compiled.plan(&structure, &policy);
        group.bench_function(BenchmarkId::new("bound_plan_sparse", label), |b| {
            b.iter(|| black_box(physical.evaluate(&structure, &policy, &groups, Some(&spatial))));
        });
    }
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_language(&mut criterion);
    bench_reuse(&mut criterion);
    criterion.final_summary();
}
