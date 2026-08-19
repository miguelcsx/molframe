//! Criterion coverage for sampled molecular-surface kernels.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, coordinates, structure};
use pdbiox_surface::{
    AtomDepthOptions, atom_depths, buried_surface, cavities, fibonacci_sphere, lee_richards,
    shrake_rupley, solvent_excluded_surface, surface_points, surface_points_at_density,
};

fn bench_sasa(c: &mut Criterion) {
    let sample_structure = structure(Sample::Small);
    let positions = coordinates(&sample_structure);
    let radii = vec![1.7_f32; positions.len()];
    let mut group = c.benchmark_group("surface_geometry");
    group.bench_function("shrake_rupley", |b| {
        b.iter(|| black_box(shrake_rupley(&positions, &radii, 1.4, 96)));
    });
    group.bench_function("lee_richards", |b| {
        b.iter(|| black_box(lee_richards(&positions, &radii, 1.4, 200)));
    });
    group.bench_function("fibonacci_sphere", |b| {
        b.iter(|| black_box(fibonacci_sphere(512)));
    });
    group.finish();

    let tiny = structure(Sample::Tiny);
    let tiny_positions = coordinates(&tiny);
    let tiny_radii = vec![1.7_f32; tiny_positions.len()];
    let membership: Vec<bool> = (0..tiny_positions.len())
        .map(|index| index < tiny_positions.len() / 2)
        .collect();
    let sampled = match surface_points(&tiny_positions, &tiny_radii, 1.4, 48) {
        Ok(points) => points,
        Err(error) => panic!("surface-point benchmark setup failed: {error}"),
    };
    let sampled_positions: Vec<[f32; 3]> = sampled.iter().map(|point| point.position).collect();
    if let Err(error) = atom_depths(
        &tiny_positions,
        &sampled_positions,
        AtomDepthOptions { cell_size: 2.5 },
    ) {
        panic!("atom-depth benchmark setup failed: {error}");
    }
    if let Err(error) = buried_surface(&tiny_positions, &tiny_radii, 1.4, 48, &membership) {
        panic!("buried-surface benchmark setup failed: {error}");
    }
    if let Err(error) = solvent_excluded_surface(&tiny_positions, &tiny_radii, 1.4, 2.0) {
        panic!("SES benchmark setup failed: {error}");
    }
    if let Err(error) = cavities(&tiny_positions, &tiny_radii, 1.4, 2.0) {
        panic!("cavity benchmark setup failed: {error}");
    }
    let mut group = c.benchmark_group("surface_extended");
    group.bench_function("surface_points", |b| {
        b.iter(|| black_box(surface_points(&tiny_positions, &tiny_radii, 1.4, 48)));
    });
    group.bench_function("surface_points_at_density", |b| {
        b.iter(|| {
            black_box(surface_points_at_density(
                &tiny_positions,
                &tiny_radii,
                1.4,
                0.3,
            ))
        });
    });
    group.bench_function("atom_depths", |b| {
        b.iter(|| {
            black_box(atom_depths(
                &tiny_positions,
                &sampled_positions,
                AtomDepthOptions { cell_size: 2.5 },
            ))
        });
    });
    group.bench_function("buried_surface", |b| {
        b.iter(|| {
            black_box(buried_surface(
                &tiny_positions,
                &tiny_radii,
                1.4,
                48,
                &membership,
            ))
        });
    });
    group.bench_function("solvent_excluded_surface", |b| {
        b.iter(|| {
            black_box(solvent_excluded_surface(
                &tiny_positions,
                &tiny_radii,
                1.4,
                2.0,
            ))
        });
    });
    group.bench_function("cavities", |b| {
        b.iter(|| black_box(cavities(&tiny_positions, &tiny_radii, 1.4, 2.0)));
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_sasa(&mut criterion);
    criterion.final_summary();
}
