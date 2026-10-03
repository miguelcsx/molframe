//! Screened contact-potential grid evaluation.

use criterion::{Criterion, black_box};
use molframe_analysis::{GridSpec, contact_potential_in};
use molframe_bench::{Sample, structure};
use molframe_core::ExecutionContext;

pub(super) fn bench_contact_potential(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let sample = structure(Sample::Small);
    let charges = vec![0.25_f64; sample.atom_count() as usize];
    let minimum = sample.positions().iter().fold([f32::MAX; 3], |low, point| {
        [
            low[0].min(point[0]),
            low[1].min(point[1]),
            low[2].min(point[2]),
        ]
    });
    let spec = GridSpec {
        voxel_to_world: [
            [1.0, 0.0, 0.0, f64::from(minimum[0]) - 8.0],
            [0.0, 1.0, 0.0, f64::from(minimum[1]) - 8.0],
            [0.0, 0.0, 1.0, f64::from(minimum[2]) - 8.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dimensions: [64, 64, 64],
    };
    let mut group = c.benchmark_group("analysis_potential");
    group.bench_function("contact_potential_64_cubed", |b| {
        b.iter(|| {
            black_box(contact_potential_in(
                &sample, &charges, spec, 12.0, &context,
            ))
        });
    });
    group.finish();
}
