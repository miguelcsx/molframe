//! Criterion coverage for spatial backends and planner dispatch.

use criterion::{Criterion, Throughput, black_box};
use pdbiox_bench::{Sample, coordinates, structure};
use pdbiox_core::execution::ExecutionContext;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::UnitCell;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, pairs_within, pairs_within_unsorted, within};

fn bench_backends(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let structure = structure(Sample::Medium);
    let positions = coordinates(&structure);
    let atom_count = match u32::try_from(positions.len()) {
        Ok(count) => count,
        Err(error) => panic!("benchmark atom count exceeds u32: {error}"),
    };
    let selection = AtomSelection::All(atom_count);
    let mut group = c.benchmark_group("spatial_pairs_within");
    group.throughput(Throughput::Elements(u64::from(atom_count)));
    for backend in [
        SpatialBackend::BruteForce,
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        group.bench_function(format!("{backend:?}"), |b| {
            b.iter(|| {
                black_box(pairs_within(
                    &positions, &selection, &selection, 4.0, backend, None, &context,
                ))
            });
        });
    }
    group.bench_function("within", |b| {
        b.iter(|| {
            black_box(within(
                &positions,
                &selection,
                &selection,
                4.0,
                SpatialBackend::Auto,
                None,
                &context,
            ))
        });
    });
    group.bench_function("same_selection_unsorted", |b| {
        b.iter(|| {
            black_box(pairs_within_unsorted(
                &positions,
                &selection,
                &selection,
                4.0,
                SpatialBackend::CellList,
                None,
                &context,
            ))
        });
    });
    group.finish();
}

fn bench_periodic_geometry(c: &mut Criterion) {
    let orthorhombic = match PeriodicBox::from_cell(UnitCell {
        lengths: [80.0, 90.0, 100.0],
        angles: [90.0, 90.0, 90.0],
    }) {
        Ok(cell) => cell,
        Err(error) => panic!("orthorhombic benchmark cell failed: {error}"),
    };
    let triclinic = match PeriodicBox::from_cell(UnitCell {
        lengths: [80.0, 90.0, 100.0],
        angles: [72.0, 81.0, 76.0],
    }) {
        Ok(cell) => cell,
        Err(error) => panic!("triclinic benchmark cell failed: {error}"),
    };
    let points: Vec<([f32; 3], [f32; 3])> = (0_u16..4_096)
        .map(|index| {
            let value = f32::from(index);
            (
                [value * 0.17, value * 0.11, value * 0.07],
                [91.0 - value * 0.13, 87.0 - value * 0.05, value * 0.19],
            )
        })
        .collect();
    let mut group = c.benchmark_group("spatial_minimum_image");
    group.throughput(Throughput::Elements(4_096));
    group.bench_function("orthorhombic", |b| {
        b.iter(|| {
            for &(left, right) in &points {
                black_box(orthorhombic.minimum_image(left, right));
            }
        });
    });
    group.bench_function("triclinic", |b| {
        b.iter(|| {
            for &(left, right) in &points {
                black_box(triclinic.minimum_image(left, right));
            }
        });
    });
    group.finish();
}

fn bench_large_scaling(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let mut group = c.benchmark_group("spatial_scaling");
    for sample in [Sample::Medium, Sample::Large] {
        let structure = structure(sample);
        let positions = coordinates(&structure);
        let selection = AtomSelection::All(structure.atom_count());
        group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
        group.bench_function(format!("cell_list/{}", sample.label()), |b| {
            b.iter(|| {
                black_box(pairs_within(
                    &positions,
                    &selection,
                    &selection,
                    4.0,
                    SpatialBackend::CellList,
                    None,
                    &context,
                ))
            });
        });
        group.bench_function(format!("auto/{}", sample.label()), |b| {
            b.iter(|| {
                black_box(pairs_within(
                    &positions,
                    &selection,
                    &selection,
                    4.0,
                    SpatialBackend::Auto,
                    None,
                    &context,
                ))
            });
        });
    }
    group.finish();
}

fn bench_streaming_reduction(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let (positions, query, target) = sparse_paired_positions(50_000);
    let mut group = c.benchmark_group("spatial_reduction");
    group.throughput(Throughput::Elements(100_000));
    group.bench_function("within_cell_list_sparse_100k", |b| {
        b.iter(|| {
            black_box(within(
                &positions,
                &query,
                &target,
                0.75,
                SpatialBackend::CellList,
                None,
                &context,
            ))
        });
    });
    group.finish();
}

fn sparse_paired_positions(pair_count: u16) -> (Vec<[f32; 3]>, AtomSelection, AtomSelection) {
    let mut positions = Vec::with_capacity(usize::from(pair_count) * 2);
    let mut query = Vec::with_capacity(usize::from(pair_count));
    let mut target = Vec::with_capacity(usize::from(pair_count));
    for pair in 0..pair_count {
        let base = f32::from(pair) * 4.0;
        let target_atom = u32::from(pair) * 2;
        positions.push([base, 0.0, 0.0]);
        positions.push([base + 0.5, 0.0, 0.0]);
        target.push(target_atom);
        query.push(target_atom + 1);
    }
    (
        positions,
        AtomSelection::from_sorted(query),
        AtomSelection::from_sorted(target),
    )
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_backends(&mut criterion);
    bench_periodic_geometry(&mut criterion);
    bench_large_scaling(&mut criterion);
    bench_streaming_reduction(&mut criterion);
    criterion.final_summary();
}
