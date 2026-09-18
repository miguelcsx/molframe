//! Exhaustive oracle and large-channel fixture for pore-profile benchmarks.

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, Throughput, black_box};
use molframe_analysis::{PoreProfileOptions, PoreSample, pore_profile};
use num_traits::ToPrimitive;
use std::f32::consts::TAU;

const MEMORY_LIMIT: usize = 100_000_000;

pub(super) fn bench_small(group: &mut BenchmarkGroup<'_, WallTime>, positions: &[[f32; 3]]) {
    let radii = vec![1.7; positions.len()];
    let options = options(16, 10.0, 1.0);
    group.bench_function("pore_profile/128", |b| {
        b.iter(|| black_box(pore_profile(positions, &radii, options)));
    });
    group.bench_function("pore_profile_exhaustive/128", |b| {
        b.iter(|| black_box(exhaustive_profile(positions, &radii, options)));
    });
}

pub(super) fn bench_scaling(c: &mut Criterion) {
    let positions = cylindrical_channel(1_000, 100);
    let radii = vec![1.7; positions.len()];
    let options = options(16, 5.0, 1.0);
    let mut group = c.benchmark_group("analysis_pore_scaling");
    group.sample_size(10);
    let element_count = match u64::try_from(positions.len()) {
        Ok(count) => count,
        Err(error) => panic!("pore benchmark size failed: {error}"),
    };
    group.throughput(Throughput::Elements(element_count));
    group.bench_function("pore_profile_exact/100000", |b| {
        b.iter(|| black_box(pore_profile(&positions, &radii, options)));
    });
    group.finish();
}

fn options(samples: usize, search_radius: f32, grid_spacing: f32) -> PoreProfileOptions {
    PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0, 0.0, 1.0],
        start: -10.0,
        end: 10.0,
        samples,
        search_radius,
        grid_spacing,
        probe_radius: 1.4,
        memory_limit_bytes: MEMORY_LIMIT,
    }
}

fn cylindrical_channel(rings: u16, atoms_per_ring: u16) -> Vec<[f32; 3]> {
    let mut positions = Vec::with_capacity(usize::from(rings) * usize::from(atoms_per_ring));
    for ring in 0..rings {
        let z = -50.0 + 100.0 * f32::from(ring) / f32::from(rings - 1);
        for atom in 0..atoms_per_ring {
            let angle = TAU * f32::from(atom) / f32::from(atoms_per_ring);
            let radius = 9.0 + 0.5 * (angle * 5.0 + z * 0.1).sin();
            positions.push([radius * angle.cos(), radius * angle.sin(), z]);
        }
    }
    positions
}

fn exhaustive_profile(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: PoreProfileOptions,
) -> Vec<PoreSample> {
    let axis = normalize(options.axis_direction);
    let first_basis = perpendicular(axis);
    let second_basis = cross(axis, first_basis);
    let Some(steps) = (2.0 * options.search_radius / options.grid_spacing)
        .ceil()
        .to_usize()
    else {
        panic!("finite pore benchmark grid")
    };
    let mut profile = Vec::with_capacity(options.samples);
    for sample in 0..options.samples {
        let fraction = if options.samples == 1 {
            0.0
        } else {
            usize_to_f32(sample) / usize_to_f32(options.samples - 1)
        };
        let axial_coordinate = options.start + fraction * (options.end - options.start);
        let axis_centre = add(options.axis_origin, scale(axis, axial_coordinate));
        let mut best_centre = axis_centre;
        let mut best_clearance = f32::NEG_INFINITY;
        for first_index in 0..=steps {
            let first = transverse_coordinate(first_index, steps, options.search_radius);
            for second_index in 0..=steps {
                let second = transverse_coordinate(second_index, steps, options.search_radius);
                if first.mul_add(first, second * second) > options.search_radius.powi(2) {
                    continue;
                }
                let candidate = add(
                    add(axis_centre, scale(first_basis, first)),
                    scale(second_basis, second),
                );
                let clearance =
                    minimum_clearance(candidate, positions, radii) - options.probe_radius;
                if clearance > best_clearance {
                    best_clearance = clearance;
                    best_centre = candidate;
                }
            }
        }
        profile.push(PoreSample {
            axial_coordinate,
            centre: best_centre,
            radius: best_clearance.max(0.0),
        });
    }
    profile
}

fn transverse_coordinate(index: usize, steps: usize, radius: f32) -> f32 {
    -radius + 2.0 * radius * usize_to_f32(index) / usize_to_f32(steps)
}

fn usize_to_f32(value: usize) -> f32 {
    match value.to_f32() {
        Some(value) => value,
        None => f32::INFINITY,
    }
}

fn minimum_clearance(point: [f32; 3], positions: &[[f32; 3]], radii: &[f32]) -> f32 {
    positions
        .iter()
        .zip(radii)
        .map(|(position, radius)| norm(subtract(point, *position)) - radius)
        .fold(f32::INFINITY, f32::min)
}

fn perpendicular(axis: [f32; 3]) -> [f32; 3] {
    let absolute = axis.map(f32::abs);
    let reference = if absolute[0] <= absolute[1] && absolute[0] <= absolute[2] {
        [1.0, 0.0, 0.0]
    } else if absolute[1] <= absolute[2] {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize(cross(axis, reference))
}

fn normalize(vector: [f32; 3]) -> [f32; 3] {
    scale(vector, norm(vector).recip())
}

fn norm(vector: [f32; 3]) -> f32 {
    vector
        .into_iter()
        .map(|value| value * value)
        .sum::<f32>()
        .sqrt()
}

fn cross(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn add(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] + right[axis])
}

fn subtract(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] - right[axis])
}

fn scale(vector: [f32; 3], factor: f32) -> [f32; 3] {
    vector.map(|value| value * factor)
}
