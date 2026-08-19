//! Criterion coverage for internal-coordinate extraction and placement.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, coordinates, structure};
use pdbiox_core::index::ModelIndex;
use pdbiox_ic::{internal_coordinates, place_atom};

fn bench_internal_coordinates(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    c.bench_function("ic_internal_coordinates/4hhb", |b| {
        b.iter(|| black_box(internal_coordinates(&structure, ModelIndex::new(0))));
    });
}

fn bench_placement(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let positions = coordinates(&structure);
    let Some(first) = positions.first().copied() else {
        panic!("tiny fixture has no first coordinate")
    };
    let Some(second) = positions.get(1).copied() else {
        panic!("tiny fixture has no second coordinate")
    };
    let Some(third) = positions.get(2).copied() else {
        panic!("tiny fixture has no third coordinate")
    };
    let points = [first, second, third];
    c.bench_function("ic_place_atom", |b| {
        b.iter(|| black_box(place_atom(points[0], points[1], points[2], 1.5, 1.9, 0.7)));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_internal_coordinates(&mut criterion);
    bench_placement(&mut criterion);
    criterion.final_summary();
}
