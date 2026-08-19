//! Criterion coverage for structural interaction and contact kernels.

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, Throughput, black_box};
use pdbiox_analysis::{
    BaseFrame, CartesianAxis, CentreGroup, DensityGridSpec, FragmentReference, GnmOptions,
    HelicalOptions, LeafletOptions, LinearDensityOptions, PoreProfileOptions,
    RadialDistributionOptions, atom_contacts, atom_contacts_between,
    atom_contacts_between_with_spatial, centre_of_mass_radial_distribution, chain_interface,
    coordination_numbers, density_map, gaussian_network_model, helical_parameters, helical_steps,
    identify_leaflets, linear_density, map_fragments, native_contact_fraction, polymer_statistics,
    pore_profile, radial_distribution, residue_contact_map, sugar_pucker, surface_contacts,
};
use pdbiox_bench::{Sample, structure};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::{Structure, contract::AnalysisPolicy};
use pdbiox_spatial::{SpatialBackend, StructureSpatial};

fn bench_contacts(c: &mut Criterion) {
    let mut group = c.benchmark_group("analysis_contacts");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium, Sample::Large] {
        let structure = structure(sample);
        group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
        group.bench_with_input(
            BenchmarkId::new("atom_contacts", sample.label()),
            &structure,
            |b, structure| {
                b.iter(|| black_box(atom_contacts(structure, 4.0, SpatialBackend::CellList)));
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
                positions, &selection, &selection, radial, None,
            ))
        });
    });
    group.bench_function("coordination_numbers", |b| {
        b.iter(|| {
            black_box(coordination_numbers(
                positions,
                &selection,
                &selection,
                2.0,
                8.0,
                SpatialBackend::CellList,
                None,
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
                positions, &selection, &selection, radial, None,
            ))
        });
    });
    group.bench_function("coordination_numbers/1aon", |b| {
        b.iter(|| {
            black_box(coordination_numbers(
                positions,
                &selection,
                &selection,
                2.0,
                8.0,
                SpatialBackend::CellList,
                None,
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
    };
    let radii = vec![1.7; positions.len()];
    let pore = PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0, 0.0, 1.0],
        start: -10.0,
        end: 10.0,
        samples: 16,
        search_radius: 10.0,
        grid_spacing: 1.0,
        probe_radius: 1.4,
    };
    let mut group = c.benchmark_group("analysis_models");
    group.bench_function("gaussian_network_model/128", |b| {
        b.iter(|| black_box(gaussian_network_model(positions, &selection, options, None)));
    });
    group.bench_function("pore_profile/128", |b| {
        b.iter(|| black_box(pore_profile(positions, &radii, pore)));
    });
    group.finish();
}

fn checked_atom_index(atom: usize) -> u32 {
    match u32::try_from(atom) {
        Ok(atom) => atom,
        Err(error) => panic!("analysis fixture atom index failed: {error}"),
    }
}

fn bench_extended_kernels(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let positions = structure.positions();
    let split = structure.atom_count() / 2;
    let left = AtomSelection::from_sorted((0..split).collect());
    let right = AtomSelection::from_sorted((split..structure.atom_count()).collect());
    let spatial = match StructureSpatial::new(
        &structure,
        &AnalysisPolicy::default(),
        SpatialBackend::CellList,
    ) {
        Ok(spatial) => spatial,
        Err(error) => panic!("analysis spatial fixture failed: {error}"),
    };
    let radii = vec![1.7; positions.len()];
    let masses = vec![12.0; positions.len()];
    let groups: Vec<_> = (0..positions.len().min(128))
        .step_by(4)
        .map(|start| CentreGroup {
            atoms: AtomSelection::from_sorted(
                (start..(start + 4).min(positions.len()))
                    .map(checked_atom_index)
                    .collect(),
            ),
        })
        .collect();
    let radial = RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 12.0,
        bins: 48,
        volume: 1_000_000.0,
        backend: SpatialBackend::CellList,
    };
    let lipid_sites = AtomSelection::from_sorted(
        (0..positions.len().min(512))
            .step_by(3)
            .map(checked_atom_index)
            .collect(),
    );
    let frames: Vec<_> = (0_u32..64)
        .map(|frame| BaseFrame {
            origin: [f64::from(frame), 0.0, 0.0],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        })
        .collect();
    let fragments = vec![
        FragmentReference {
            id: "fragment-a".into(),
            coordinates: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        },
        FragmentReference {
            id: "fragment-b".into(),
            coordinates: vec![
                [0.0, 0.0, 0.0],
                [1.1, 0.0, 0.0],
                [0.0, 1.1, 0.0],
                [0.0, 0.0, 1.1],
            ],
        },
    ];
    let trace: Vec<_> = (0_u16..64)
        .map(|index| {
            Some([
                f32::from(index % 4),
                f32::from((index + 1) % 4),
                f32::from((index + 2) % 4),
            ])
        })
        .collect();
    let mut group = c.benchmark_group("analysis_extended");
    bench_extended_contacts(&mut group, &structure, &left, &right, &spatial);
    bench_extended_surface(&mut group, &structure, &radii);
    bench_extended_distributions(
        &mut group,
        positions,
        &masses,
        &groups,
        radial,
        &lipid_sites,
    );
    bench_extended_reference_kernels(&mut group, &frames, &trace, &fragments);
    group.finish();
}

fn bench_extended_contacts(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &Structure,
    left: &AtomSelection,
    right: &AtomSelection,
    spatial: &StructureSpatial<'_>,
) {
    group.bench_function("atom_contacts_between/half", |b| {
        b.iter(|| {
            black_box(atom_contacts_between(
                structure,
                left,
                right,
                4.0,
                SpatialBackend::CellList,
            ))
        });
    });
    group.bench_function("atom_contacts_between_with_spatial/half", |b| {
        b.iter(|| {
            black_box(atom_contacts_between_with_spatial(
                structure,
                left,
                right,
                4.0,
                SpatialBackend::CellList,
                spatial,
            ))
        });
    });
}

fn bench_extended_surface(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &Structure,
    radii: &[f32],
) {
    group.bench_function("surface_contacts/medium", |b| {
        b.iter(|| {
            black_box(surface_contacts(
                structure,
                radii,
                0.4,
                1.4,
                1.0,
                0.5,
                SpatialBackend::CellList,
            ))
        });
    });
}

fn bench_extended_distributions(
    group: &mut BenchmarkGroup<'_, WallTime>,
    positions: &[[f32; 3]],
    masses: &[f64],
    groups: &[CentreGroup],
    radial: RadialDistributionOptions,
    lipid_sites: &AtomSelection,
) {
    group.bench_function("centre_of_mass_radial/32", |b| {
        b.iter(|| {
            black_box(centre_of_mass_radial_distribution(
                positions, masses, groups, groups, radial, None,
            ))
        });
    });
    group.bench_function("identify_leaflets/170", |b| {
        b.iter(|| {
            black_box(identify_leaflets(
                positions,
                lipid_sites,
                LeafletOptions {
                    connection_distance: 6.0,
                    backend: SpatialBackend::CellList,
                },
                None,
            ))
        });
    });
}

fn bench_extended_reference_kernels(
    group: &mut BenchmarkGroup<'_, WallTime>,
    frames: &[BaseFrame],
    trace: &[Option<[f32; 3]>],
    fragments: &[FragmentReference],
) {
    group.bench_function("map_fragments/64x4", |b| {
        b.iter(|| black_box(map_fragments(trace, fragments, 2.0)));
    });
    group.bench_function("helical_parameters", |b| {
        b.iter(|| {
            black_box(helical_parameters(
                frames[0],
                frames[1],
                HelicalOptions {
                    frame_tolerance: 1.0e-12,
                },
            ))
        });
    });
    group.bench_function("helical_steps/64", |b| {
        b.iter(|| {
            black_box(helical_steps(
                frames,
                HelicalOptions {
                    frame_tolerance: 1.0e-12,
                },
            ))
        });
    });
    group.bench_function("sugar_pucker", |b| {
        b.iter(|| black_box(sugar_pucker([12.0, -8.0, 5.0, 16.0, -4.0])));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_contacts(&mut criterion);
    bench_physical_kernels(&mut criterion);
    bench_large_physical_scaling(&mut criterion);
    bench_model_kernels(&mut criterion);
    bench_extended_kernels(&mut criterion);
    criterion.final_summary();
}
