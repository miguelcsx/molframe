//! Extended analysis benchmark families.

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, black_box};
use molframe_analysis::{
    BaseFrame, CentreGroup, FragmentReference, HelicalOptions, LeafletOptions,
    RadialDistributionOptions, atom_contacts_between, atom_contacts_between_with_spatial,
    centre_of_mass_radial_distribution, helical_parameters, helical_steps, identify_leaflets,
    map_fragments, sugar_pucker, surface_contacts,
};
use molframe_bench::{Sample, structure};
use molframe_core::selection::AtomSelection;
use molframe_core::{ExecutionContext, Structure, contract::AnalysisPolicy};
use molframe_spatial::{SpatialBackend, StructureSpatial};

fn context() -> ExecutionContext {
    ExecutionContext::default()
}

fn checked_atom_index(atom: usize) -> u32 {
    match u32::try_from(atom) {
        Ok(atom) => atom,
        Err(error) => panic!("analysis fixture atom index failed: {error}"),
    }
}

pub(super) fn bench_extended_kernels(c: &mut Criterion) {
    let context = context();
    let structure = structure(Sample::Medium);
    let positions = structure.positions();
    let split = structure.atom_count() / 2;
    let left = AtomSelection::from_sorted((0..split).collect());
    let right = AtomSelection::from_sorted((split..structure.atom_count()).collect());
    let spatial = match StructureSpatial::new(
        &structure,
        &AnalysisPolicy::default(),
        SpatialBackend::CellList,
        &context,
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
    let fragments = fragment_references();
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

fn fragment_references() -> Vec<FragmentReference> {
    vec![
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
    ]
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
                &context(),
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
                molframe_analysis::SurfaceContactOptions {
                    tolerance: 0.4,
                    probe: 1.4,
                    surface_density: 1.0,
                    minimum_area: 0.5,
                    backend: SpatialBackend::CellList,
                },
                &context(),
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
                positions,
                masses,
                groups,
                groups,
                radial,
                None,
                &context(),
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
                &context(),
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
