//! Criterion coverage for structural interaction and contact kernels.

use criterion::{BenchmarkId, Criterion, Throughput, black_box};
use molframe_analysis::{
    CartesianAxis, DensityGridSpec, GnmOptions, LinearDensityOptions, RadialDistributionOptions,
    atom_contacts, chain_interface, coordination_numbers, density_map, gaussian_network_model,
    linear_density, native_contact_fraction, polymer_statistics, radial_distribution,
    residue_contact_map,
};
use molframe_bench::{Sample, structure};
use molframe_core::ExecutionContext;
use molframe_core::selection::AtomSelection;
use molframe_spatial::SpatialBackend;

#[path = "analysis/extended.rs"]
mod extended;
use extended::{bench_extended_kernels, bench_hydrogen_bonds};

#[path = "analysis/gnm.rs"]
mod gnm;
use gnm::{dense_gnm_reference, grid as gnm_grid};
#[path = "analysis/pore.rs"]
mod pore;

fn context() -> ExecutionContext {
    ExecutionContext::default()
}

fn bench_contacts(c: &mut Criterion) {
    let mut group = c.benchmark_group("analysis_contacts");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium, Sample::Large] {
        let structure = structure(sample);
        group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
        group.bench_with_input(
            BenchmarkId::new("atom_contacts", sample.label()),
            &structure,
            |b, structure| {
                b.iter(|| {
                    black_box(atom_contacts(
                        structure,
                        4.0,
                        SpatialBackend::CellList,
                        &context(),
                    ))
                });
            },
        );
    }

    let structure = structure(Sample::Medium);
    group.bench_function("residue_contact_map", |b| {
        b.iter(|| {
            black_box(residue_contact_map(
                &structure,
                6.0,
                2,
                SpatialBackend::CellList,
                &context(),
            ))
        });
    });
    group.bench_function("native_contact_fraction", |b| {
        b.iter(|| {
            black_box(native_contact_fraction(
                &structure,
                &structure,
                4.0,
                1.2,
                SpatialBackend::CellList,
                &context(),
            ))
        });
    });
    group.bench_function("chain_interface/A-B", |b| {
        b.iter(|| {
            black_box(chain_interface(
                &structure,
                "A",
                "B",
                5.0,
                SpatialBackend::CellList,
                &context(),
            ))
        });
    });
    group.finish();
}

fn bench_physical_kernels(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let positions = structure.positions();
    let selection = AtomSelection::All(structure.atom_count());
    let weights = vec![1.0; positions.len()];
    let radial = RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 12.0,
        bins: 48,
        volume: 1_000_000.0,
        backend: SpatialBackend::CellList,
    };
    let mut group = c.benchmark_group("analysis_physical");
    group.bench_function("radial_distribution", |b| {
        b.iter(|| {
            black_box(radial_distribution(
                positions,
                &selection,
                &selection,
                radial,
                None,
                &context(),
            ))
        });
    });
    group.bench_function("coordination_numbers", |b| {
        b.iter(|| {
            black_box(coordination_numbers(
                positions,
                &selection,
                &selection,
                molframe_analysis::CoordinationOptions {
                    minimum_distance: 2.0,
                    maximum_distance: 8.0,
                    backend: SpatialBackend::CellList,
                },
                None,
                &context(),
            ))
        });
    });
    group.bench_function("linear_density", |b| {
        b.iter(|| {
            black_box(linear_density(
                positions,
                &weights,
                LinearDensityOptions {
                    axis: CartesianAxis::Z,
                    minimum: -200.0,
                    maximum: 200.0,
                    bins: 128,
                },
            ))
        });
    });
    group.bench_function("density_map", |b| {
        b.iter(|| {
            black_box(density_map(
                positions,
                &weights,
                DensityGridSpec {
                    origin: [-200.0; 3],
                    spacing: [5.0; 3],
                    shape: [80; 3],
                },
            ))
        });
    });
    group.bench_function("polymer_statistics", |b| {
        b.iter(|| black_box(polymer_statistics(&positions[..positions.len().min(1_024)])));
    });
    group.finish();
}

fn bench_large_physical_scaling(c: &mut Criterion) {
    let structure = structure(Sample::Large);
    let positions = structure.positions();
    let selection = AtomSelection::All(structure.atom_count());
    let weights = vec![1.0; positions.len()];
    let radial = RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 12.0,
        bins: 48,
        volume: 1_000_000.0,
        backend: SpatialBackend::CellList,
    };
    let mut group = c.benchmark_group("analysis_scaling");
    group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
    group.bench_function("radial_distribution/1aon", |b| {
        b.iter(|| {
            black_box(radial_distribution(
                positions,
                &selection,
                &selection,
                radial,
                None,
                &context(),
            ))
        });
    });
    group.bench_function("coordination_numbers/1aon", |b| {
        b.iter(|| {
            black_box(coordination_numbers(
                positions,
                &selection,
                &selection,
                molframe_analysis::CoordinationOptions {
                    minimum_distance: 2.0,
                    maximum_distance: 8.0,
                    backend: SpatialBackend::CellList,
                },
                None,
                &context(),
            ))
        });
    });
    group.bench_function("linear_density/1aon", |b| {
        b.iter(|| {
            black_box(linear_density(
                positions,
                &weights,
                LinearDensityOptions {
                    axis: CartesianAxis::Z,
                    minimum: -200.0,
                    maximum: 200.0,
                    bins: 128,
                },
            ))
        });
    });
    group.bench_function("density_map/1aon", |b| {
        b.iter(|| {
            black_box(density_map(
                positions,
                &weights,
                DensityGridSpec {
                    origin: [-200.0; 3],
                    spacing: [5.0; 3],
                    shape: [80; 3],
                },
            ))
        });
    });
    group.finish();
}

fn bench_model_kernels(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let count = structure.positions().len().min(128);
    let positions = &structure.positions()[..count];
    let count_u32 = match u32::try_from(count) {
        Ok(count) => count,
        Err(error) => panic!("analysis model benchmark size failed: {error}"),
    };
    let selection = AtomSelection::from_sorted((0..count_u32).collect());
    let options = GnmOptions {
        contact_distance: 10.0,
        mode_count: 8,
        zero_mode_tolerance: 1e-8,
        memory_limit_bytes: 16 * 1_024 * 1_024,
        backend: SpatialBackend::CellList,
        reduction: molframe_core::parallel::ReductionPolicy::Deterministic,
    };
    let mut group = c.benchmark_group("analysis_models");
    group.bench_function("gaussian_network_model/128", |b| {
        b.iter(|| {
            black_box(gaussian_network_model(
                positions,
                &selection,
                options,
                None,
                &context(),
            ))
        });
    });
    group.bench_function("gaussian_network_model_dense_reference/128", |b| {
        b.iter(|| black_box(dense_gnm_reference(positions, &selection, options)));
    });
    pore::bench_small(&mut group, positions);
    group.finish();

    let grid = gnm_grid(25, 20, 20);
    let grid_selection = AtomSelection::All(10_000);
    let grid_options = GnmOptions {
        contact_distance: 1.01,
        mode_count: 1,
        zero_mode_tolerance: 1e-10,
        memory_limit_bytes: 500_000_000,
        backend: SpatialBackend::CellList,
        reduction: molframe_core::parallel::ReductionPolicy::Deterministic,
    };
    let mut scaling = c.benchmark_group("analysis_gnm_scaling");
    scaling.sample_size(10);
    scaling.throughput(Throughput::Elements(10_000));
    scaling.bench_function("gaussian_network_model_sparse/10000", |b| {
        b.iter(|| {
            black_box(gaussian_network_model(
                &grid,
                &grid_selection,
                grid_options,
                None,
                &context(),
            ))
        });
    });
    scaling.finish();
    pore::bench_scaling(c);
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_contacts(&mut criterion);
    bench_physical_kernels(&mut criterion);
    bench_large_physical_scaling(&mut criterion);
    bench_model_kernels(&mut criterion);
    bench_extended_kernels(&mut criterion);
    bench_hydrogen_bonds(&mut criterion);
    criterion.final_summary();
}
