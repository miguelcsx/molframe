//! Criterion coverage for spatial backends and planner dispatch.

use criterion::{Criterion, Throughput, black_box};
use molframe_bench::{Sample, coordinates, structure};
use molframe_core::execution::ExecutionContext;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::UnitCell;
use molframe_spatial::{PeriodicBox, SpatialBackend, pairs_within, pairs_within_unsorted, within};

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

/// Uniform points filling `lengths`, drawn from the shared `Seed` sequence so
/// every run sees the same fixture.
fn periodic_positions(lengths: [f64; 3], atom_count: u32) -> Vec<[f32; 3]> {
    use num_traits::ToPrimitive;

    fn f32_lossless(value: f64) -> f32 {
        let Some(converted) = value.to_f32() else {
            return f32::INFINITY;
        };
        converted
    }

    let mut seed = molframe_bench::Seed::new(0x5EED);
    (0..atom_count)
        .map(|_| {
            [
                f32_lossless(seed.next_unit() * lengths[0]),
                f32_lossless(seed.next_unit() * lengths[1]),
                f32_lossless(seed.next_unit() * lengths[2]),
            ]
        })
        .collect()
}

fn periodic_box(angles: [f64; 3], lengths: [f64; 3]) -> PeriodicBox {
    match PeriodicBox::from_cell(UnitCell { lengths, angles }) {
        Ok(box_) => box_,
        Err(error) => panic!("periodic benchmark cell failed: {error}"),
    }
}

/// Periodic `pairs_within` at the 19-performance budget scale (100k atoms,
/// 5 Å): today the planner pins this to scalar brute force, so the per-backend
/// rows are the before/after evidence for the cost-model fix.
fn bench_periodic_search(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let orthorhombic = periodic_box([90.0; 3], [120.0; 3]);
    let triclinic = periodic_box([72.0, 81.0, 76.0], [120.0; 3]);
    let mut group = c.benchmark_group("spatial_pairs_within_periodic");
    for (atom_count, label) in [(10_000_u32, "10k"), (100_000, "100k")] {
        let positions = periodic_positions([120.0; 3], atom_count);
        let selection = AtomSelection::All(atom_count);
        group.throughput(Throughput::Elements(u64::from(atom_count)));
        for (cell_name, cell) in [("orthorhombic", &orthorhombic), ("triclinic", &triclinic)] {
            for backend in [
                SpatialBackend::BruteForce,
                SpatialBackend::CellList,
                SpatialBackend::KdTree,
            ] {
                group.bench_function(format!("{label}/{cell_name}/{backend:?}"), |b| {
                    b.iter(|| {
                        black_box(pairs_within(
                            &positions,
                            &selection,
                            &selection,
                            5.0,
                            backend,
                            Some(cell),
                            &context,
                        ))
                    });
                });
            }
        }
    }
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_backends(&mut criterion);
    bench_periodic_geometry(&mut criterion);
    bench_periodic_search(&mut criterion);
    bench_large_scaling(&mut criterion);
    bench_streaming_reduction(&mut criterion);
    criterion.final_summary();
}
