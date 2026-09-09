//! Criterion coverage for structure-quality validation kernels.

use criterion::{BenchmarkGroup, Criterion, black_box, measurement::WallTime};
use pdbiox_bench::{Sample, structure};
use pdbiox_chem::RadiusSet;
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::Namespace;
use pdbiox_core::index::AtomIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::UnitCell;
use pdbiox_geom::EigenOptions;
use pdbiox_spatial::SpatialBackend;
use pdbiox_validate::{
    AltlocOccupancyOptions, BondDeviation, PlanarityOptions, PlaneRestraint, ReferenceDistribution,
    ReferenceLibrary, TlsGroup, TlsModel, altloc_occupancy_sums, assess_bond_deviation,
    b_factor_distribution, bond_length_deviations, chirality_outliers, cis_peptides, clashes,
    completeness, ligand_geometry, ligand_geometry_outliers, masked_real_space_correlation,
    nonplanar_aromatic_rings, overvalent_atoms, plane_restraint_outliers, quality_flags,
    real_space_map_correlation, sampled_real_space_correlation, tls_b_factor_consistency,
};
use pdbiox_xtal::{DensityMap, MapBoundary};

fn bench_validation(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let mut group = c.benchmark_group("validate_structure");
    group.bench_function("clashes", |b| {
        b.iter(|| {
            black_box(clashes(
                &structure,
                0.4,
                RadiusSet::Bondi,
                SpatialBackend::CellList,
                &ExecutionContext::default(),
            ))
        });
    });
    group.bench_function("occupancy", |b| {
        let options = AltlocOccupancyOptions {
            expected_sum: 1.0,
            tolerance: 0.02,
        };
        b.iter(|| black_box(altloc_occupancy_sums(&structure, Namespace::Auth, options)));
    });
    group.bench_function("completeness", |b| {
        b.iter(|| black_box(completeness(&structure, Namespace::Auth)));
    });
    let selection = AtomSelection::All(structure.atom_count());
    group.bench_function("b_factor_distribution", |b| {
        b.iter(|| black_box(b_factor_distribution(&structure, &selection, 3.0)));
    });
    group.finish();
}

fn bench_extended_validation(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let mut group = c.benchmark_group("validate_extended");
    bench_extended_structure(&mut group, &structure);
    bench_extended_geometry(&mut group, &structure);
    bench_extended_maps(&mut group, &structure);
    bench_extended_references(&mut group);
    bench_extended_chirality(&mut group, &structure);
    group.finish();
}

fn bench_extended_structure(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &pdbiox_core::structure::Structure,
) {
    group.bench_function("quality_flags/medium", |b| {
        b.iter(|| black_box(quality_flags(structure)));
    });
    group.bench_function("bond_length_deviations/medium", |b| {
        b.iter(|| black_box(bond_length_deviations(structure, 0.2)));
    });
    group.bench_function("overvalent_atoms/medium", |b| {
        b.iter(|| black_box(overvalent_atoms(structure)));
    });
    group.bench_function("cis_peptides/medium", |b| {
        b.iter(|| black_box(cis_peptides(structure, 30.0)));
    });
    group.bench_function("ligand_geometry/medium", |b| {
        b.iter(|| black_box(ligand_geometry(structure, 0.2)));
    });
    group.bench_function("ligand_geometry_outliers/medium", |b| {
        b.iter(|| black_box(ligand_geometry_outliers(structure, 0.2)));
    });
}

fn bench_extended_geometry(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &pdbiox_core::structure::Structure,
) {
    let plane_selection =
        AtomSelection::from_sorted((0..structure.atom_count().min(128)).collect());
    let plane_options = PlanarityOptions {
        maximum_deviation: 0.05,
        plane_fit: EigenOptions::standard(),
    };
    group.bench_function("nonplanar_aromatic_rings/medium", |b| {
        b.iter(|| black_box(nonplanar_aromatic_rings(structure, plane_options)));
    });
    group.bench_function("plane_restraint_outliers/128", |b| {
        let restraints = [PlaneRestraint {
            id: "all".into(),
            atoms: plane_selection.clone(),
        }];
        b.iter(|| {
            black_box(plane_restraint_outliers(
                structure,
                &restraints,
                plane_options,
            ))
        });
    });
    group.bench_function("tls_b_factor_consistency/128", |b| {
        let tls_group = TlsGroup {
            id: "all".into(),
            atoms: plane_selection.clone(),
            model: TlsModel {
                origin: [0.0; 3],
                translation: [[0.0; 3]; 3],
                libration: [[0.0; 3]; 3],
                screw: [[0.0; 3]; 3],
            },
        };
        b.iter(|| {
            black_box(tls_b_factor_consistency(
                structure,
                std::slice::from_ref(&tls_group),
                1.0,
                1.0e-12,
            ))
        });
    });
}

fn bench_extended_maps(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &pdbiox_core::structure::Structure,
) {
    let map = density_map();
    let mask: Vec<_> = (0..map.values.len()).map(|index| index % 3 != 0).collect();
    let sample_positions: Vec<[f64; 3]> = structure
        .positions()
        .iter()
        .take(512)
        .map(|position| position.map(f64::from))
        .collect();
    group.bench_function("real_space_map_correlation/4096", |b| {
        b.iter(|| black_box(real_space_map_correlation(&map, &map)));
    });
    group.bench_function("masked_real_space_map_correlation/4096", |b| {
        b.iter(|| black_box(masked_real_space_correlation(&map, &map, &mask)));
    });
    group.bench_function("sampled_real_space_map_correlation/512", |b| {
        b.iter(|| {
            black_box(sampled_real_space_correlation(
                &map,
                &map,
                &sample_positions,
                MapBoundary::Missing,
            ))
        });
    });
}

fn bench_extended_references(group: &mut BenchmarkGroup<'_, WallTime>) {
    let distribution = match ReferenceDistribution::histogram(
        "bond",
        vec![-2.0, -1.0, 0.0, 1.0, 2.0],
        vec![1.0, 2.0, 4.0, 2.0],
    ) {
        Ok(distribution) => distribution,
        Err(error) => panic!("validation reference fixture failed: {error}"),
    };
    let references = match ReferenceLibrary::new("bench", "v1", [distribution]) {
        Ok(references) => references,
        Err(error) => panic!("validation reference library failed: {error}"),
    };
    let deviation = BondDeviation {
        atom_a: AtomIndex::new(0),
        atom_b: AtomIndex::new(1),
        observed: 1.5,
        expected: 1.0,
        deviation: 0.5,
    };
    group.bench_function("reference_assess_bond", |b| {
        b.iter(|| black_box(assess_bond_deviation(&deviation, &references, "bond")));
    });
}

fn bench_extended_chirality(
    group: &mut BenchmarkGroup<'_, WallTime>,
    structure: &pdbiox_core::structure::Structure,
) {
    let provider = match pdbiox_chem::MemoryProvider::new(
        pdbiox_core::contract::DictionaryVersion::new("empty"),
        std::iter::empty::<pdbiox_chem::Component>(),
    ) {
        Ok(provider) => provider,
        Err(error) => panic!("empty chemistry provider failed: {error}"),
    };
    group.bench_function("chirality_outliers/medium", |b| {
        b.iter(|| {
            black_box(chirality_outliers(
                structure,
                &provider,
                &pdbiox_core::contract::AnalysisPolicy::default(),
                pdbiox_validate::ChiralityOptions {
                    minimum_abs_volume: 1.0e-6,
                },
            ))
        });
    });
}

fn density_map() -> DensityMap {
    let values: Vec<f32> = (0_u16..4_096)
        .map(|index| f32::from(index % 97) * 0.1)
        .collect();
    DensityMap {
        dimensions: [16, 16, 16],
        starts: [0, 0, 0],
        sampling: [16, 16, 16],
        cell: UnitCell {
            lengths: [16.0, 16.0, 16.0],
            angles: [90.0, 90.0, 90.0],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values,
    }
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_validation(&mut criterion);
    bench_extended_validation(&mut criterion);
    criterion.final_summary();
}
