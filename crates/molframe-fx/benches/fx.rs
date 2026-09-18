//! Criterion coverage for declarative functional-geometry evaluation.

use criterion::{Criterion, black_box};
use molframe_bench::{Sample, structure};
use molframe_fx::{ame_heavy_atom_1_0, motifbench_1_0};
use std::collections::BTreeMap;

fn bench_profiles(c: &mut Criterion) {
    let metrics = BTreeMap::from([
        (Box::<str>::from("rmsd"), 0.25),
        (Box::<str>::from("tm_score"), 0.8),
    ]);
    let mut group = c.benchmark_group("fx_profile_decisions");
    group.bench_function("motifbench", |b| {
        b.iter(|| black_box(motifbench_1_0().decide_candidate(&metrics)));
    });
    group.bench_function("ame", |b| {
        b.iter(|| black_box(ame_heavy_atom_1_0().decide_candidate(&metrics)));
    });
    group.finish();
}

fn bench_structure_input(c: &mut Criterion) {
    let structure = structure(Sample::Small);
    c.bench_function("fx_structure_snapshot", |b| {
        b.iter(|| black_box(structure.clone()));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_profiles(&mut criterion);
    bench_structure_input(&mut criterion);
    criterion.final_summary();
}
