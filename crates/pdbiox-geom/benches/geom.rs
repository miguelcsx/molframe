//! Criterion coverage for coordinate and matrix geometry kernels.

use criterion::{Criterion, Throughput, black_box};
use pdbiox_bench::{Sample, coordinates, perturbed, structure};
use pdbiox_geom::{
    angle, asphericity, best_fit_plane, centroid, dihedral, distance_matrix, inertia_tensor,
    principal_axes, radius_of_gyration, rmsd, rmsf, superpose,
};

fn bench_coordinate_kernels(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let reference = coordinates(&structure);
    let model = perturbed(&reference, 0.02);
    let frames = [&reference[..], &model[..]];
    let mut group = c.benchmark_group("geom_coordinates");
    group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
    group.bench_function("centroid", |b| b.iter(|| black_box(centroid(&reference))));
    group.bench_function("distance_matrix", |b| {
        b.iter(|| black_box(distance_matrix(&reference)));
    });
    group.bench_function("rmsd", |b| b.iter(|| black_box(rmsd(&model, &reference))));
    group.bench_function("superpose", |b| {
        b.iter(|| black_box(superpose(&model, &reference)));
    });
    group.bench_function("rmsf", |b| b.iter(|| black_box(rmsf(&frames))));
    group.finish();

    let masses = vec![12.0; reference.len()];
    let mut group = c.benchmark_group("geom_moments");
    group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
    group.bench_function("radius_of_gyration", |b| {
        b.iter(|| black_box(radius_of_gyration(&reference, &masses)));
    });
    group.bench_function("inertia_tensor", |b| {
        b.iter(|| black_box(inertia_tensor(&reference, &masses)));
    });
    group.bench_function("principal_axes", |b| {
        b.iter(|| black_box(principal_axes(&reference, &masses)));
    });
    group.bench_function("asphericity", |b| {
        b.iter(|| black_box(asphericity(&reference)));
    });
    group.bench_function("best_fit_plane", |b| {
        b.iter(|| black_box(best_fit_plane(&reference)));
    });
    group.finish();

    let quadruples: Vec<_> = reference.windows(4).collect();
    c.bench_function("geom_local/angles_and_dihedrals", |b| {
        b.iter(|| {
            for points in &quadruples {
                black_box(angle(points[0], points[1], points[2]));
                black_box(dihedral(points[0], points[1], points[2], points[3]));
            }
        });
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_coordinate_kernels(&mut criterion);
    criterion.final_summary();
}
