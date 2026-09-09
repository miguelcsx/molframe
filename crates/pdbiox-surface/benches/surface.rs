//! Criterion coverage for sampled molecular-surface kernels.

use criterion::{Criterion, Throughput, black_box};
use pdbiox_bench::{Sample, coordinates, structure};
use pdbiox_core::ExecutionContext;
use pdbiox_surface::{
    AtomDepthOptions, SurfaceGridOptions, atom_depths, buried_surface, cavities,
    cavities_with_options, fibonacci_sphere, lee_richards, shrake_rupley, solvent_excluded_surface,
    solvent_excluded_surface_with_options, surface_points, surface_points_at_density,
};

fn bench_sasa(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let sample_structure = structure(Sample::Small);
    let positions = coordinates(&sample_structure);
    let radii = vec![1.7_f32; positions.len()];
    let mut group = c.benchmark_group("surface_geometry");
    group.bench_function("shrake_rupley", |b| {
        b.iter(|| black_box(shrake_rupley(&positions, &radii, 1.4, 96, &context)));
    });
    group.bench_function("lee_richards", |b| {
        b.iter(|| black_box(lee_richards(&positions, &radii, 1.4, 200, &context)));
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
    let sampled = match surface_points(&tiny_positions, &tiny_radii, 1.4, 48, &context) {
        Ok(points) => points,
        Err(error) => panic!("surface-point benchmark setup failed: {error}"),
    };
    let sampled_positions: Vec<[f32; 3]> = sampled.iter().map(|point| point.position).collect();
    validate_extended_setup(
        &tiny_positions,
        &tiny_radii,
        &sampled_positions,
        &membership,
        &context,
    );
    let mut group = c.benchmark_group("surface_extended");
    group.bench_function("surface_points", |b| {
        b.iter(|| {
            black_box(surface_points(
                &tiny_positions,
                &tiny_radii,
                1.4,
                48,
                &context,
            ))
        });
    });
    group.bench_function("surface_points_at_density", |b| {
        b.iter(|| {
            black_box(surface_points_at_density(
                &tiny_positions,
                &tiny_radii,
                1.4,
                0.3,
                &context,
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
                &context,
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

fn validate_extended_setup(
    positions: &[[f32; 3]],
    radii: &[f32],
    sampled: &[[f32; 3]],
    membership: &[bool],
    context: &ExecutionContext,
) {
    if let Err(error) = atom_depths(positions, sampled, AtomDepthOptions { cell_size: 2.5 }) {
        panic!("atom-depth benchmark setup failed: {error}");
    }
    if let Err(error) = buried_surface(positions, radii, 1.4, 48, membership, context) {
        panic!("buried-surface benchmark setup failed: {error}");
    }
    if let Err(error) = solvent_excluded_surface(positions, radii, 1.4, 2.0) {
        panic!("SES benchmark setup failed: {error}");
    }
    if let Err(error) = cavities(positions, radii, 1.4, 2.0) {
        panic!("cavity benchmark setup failed: {error}");
    }
}

fn bench_voxel_stress(c: &mut Criterion) {
    // A single sphere of radius 4.85 A at 0.1 A resolution produces exactly a
    // 100 x 100 x 100 grid. This isolates the million-cell flood, distance
    // transform and polygonisation paths from structure-parsing overhead.
    let stress_positions = [[0.0_f32; 3]];
    let stress_radii = [4.85_f32];
    let stress_options = SurfaceGridOptions {
        resolution: 0.1,
        max_cells: 1_000_000,
        max_workspace_bytes: SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
    };
    if let Err(error) = cavities_with_options(&stress_positions, &stress_radii, 0.0, stress_options)
    {
        panic!("million-cell cavity benchmark setup failed: {error}");
    }
    if let Err(error) =
        solvent_excluded_surface_with_options(&stress_positions, &stress_radii, 0.0, stress_options)
    {
        panic!("million-cell SES benchmark setup failed: {error}");
    }
    let mut group = c.benchmark_group("surface_voxel_stress_1m_cells");
    group.sample_size(10);
    group.throughput(Throughput::Elements(1_000_000));
    group.bench_function("solvent_excluded_surface", |b| {
        b.iter(|| {
            black_box(solvent_excluded_surface_with_options(
                &stress_positions,
                &stress_radii,
                0.0,
                stress_options,
            ))
        });
    });
    group.bench_function("cavities", |b| {
        b.iter(|| {
            black_box(cavities_with_options(
                &stress_positions,
                &stress_radii,
                0.0,
                stress_options,
            ))
        });
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_sasa(&mut criterion);
    bench_voxel_stress(&mut criterion);
    criterion.final_summary();
}
