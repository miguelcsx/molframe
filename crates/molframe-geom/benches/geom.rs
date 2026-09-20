//! Criterion coverage for coordinate and matrix geometry kernels.

use criterion::{Criterion, Throughput, black_box};
use molframe_bench::{Sample, coordinates, perturbed, structure};
use molframe_geom::{
    angle, angles_into, asphericity, best_fit_plane, centroid, dihedral, distance, distance_matrix,
    distance_matrix_into, distances_into, inertia_tensor, principal_axes, radius_of_gyration, rmsd,
    rmsf, superpose, torsions_into,
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
    let mut matrix = vec![0.0; reference.len() * reference.len()];
    group.bench_function("distance_matrix_into", |b| {
        b.iter(|| black_box(distance_matrix_into(&reference, &mut matrix)));
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

    bench_batch_measurements(c, &reference, &model);
}

fn bench_batch_measurements(c: &mut Criterion, reference: &[[f32; 3]], model: &[[f32; 3]]) {
    let rows = reference.len().saturating_sub(3);
    let first = &reference[..rows];
    let second = &reference[1..=rows];
    let third = &reference[2..=rows + 1];
    let fourth = &reference[3..=rows + 2];
    let paired = &model[..rows];
    let mut output = vec![0.0; rows];
    let mut group = c.benchmark_group("geom_batch");
    let rows_u64 = match u64::try_from(rows) {
        Ok(rows) => rows,
        Err(_) => u64::MAX,
    };
    group.throughput(Throughput::Elements(rows_u64));
    group.bench_function("distances_scalar", |b| {
        b.iter(|| {
            for row in 0..rows {
                output[row] = distance(first[row], paired[row]);
            }
            black_box(&output);
        });
    });
    group.bench_function("distances_simd_into", |b| {
        b.iter(|| black_box(distances_into(first, paired, &mut output)));
    });
    group.bench_function("angles_scalar", |b| {
        b.iter(|| {
            for row in 0..rows {
                output[row] = match angle(first[row], second[row], third[row]) {
                    Some(angle) => angle,
                    None => f64::NAN,
                };
            }
            black_box(&output);
        });
    });
    group.bench_function("angles_simd_into", |b| {
        b.iter(|| black_box(angles_into(first, second, third, &mut output)));
    });
    group.bench_function("torsions_scalar", |b| {
        b.iter(|| {
            for row in 0..rows {
                output[row] = match dihedral(first[row], second[row], third[row], fourth[row]) {
                    Some(dihedral) => dihedral,
                    None => f64::NAN,
                };
            }
            black_box(&output);
        });
    });
    group.bench_function("torsions_simd_into", |b| {
        b.iter(|| black_box(torsions_into(first, second, third, fourth, &mut output)));
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_coordinate_kernels(&mut criterion);
    criterion.final_summary();
}
