//! Extended analysis benchmark families.

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, Throughput, black_box};
use molframe_analysis::{
    BaseFrame, CentreGroup, FragmentReference, HelicalOptions, HydrogenBondOptions, LeafletOptions,
    RadialDistributionOptions, atom_contacts_between, atom_contacts_between_with_spatial,
    centre_of_mass_radial_distribution, helical_parameters, helical_steps, hydrogen_bonds,
    identify_leaflets, map_fragments, sugar_pucker, surface_contacts,
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

/// A water lattice of `side^3` O–H pairs: every oxygen both donor and
/// acceptor, spacing 2.8 ångström. The same shape the hbond property tests
/// use, scaled to the kernel's block boundaries. Even indices are oxygens.
fn hydrogen_bond_lattice(side: i16) -> Structure {
    use molframe_core::io::{InputBuffer, ReadOptions};
    use molframe_core::{
        AnnotationColumn, AtomAnnotation, BondOrder, BondProvenance, BondRecord, BondTableBuilder,
        Presence,
    };
    use std::fmt::Write;

    let mut source = String::from(
        "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut serial = 0_u32;
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                let (ox, oy, oz) = (f32::from(x) * 2.8, f32::from(y) * 2.8, f32::from(z) * 2.8);
                let residue = serial / 2 + 1;
                serial += 1;
                writeln!(
                    source,
                    "ATOM {serial} O O HOH A {residue} {ox:.3} {oy:.3} {oz:.3}"
                )
                .expect("fixture row");
                serial += 1;
                writeln!(
                    source,
                    "ATOM {serial} H H1 HOH A {residue} {:.3} {oy:.3} {oz:.3}",
                    ox + 0.96
                )
                .expect("fixture row");
            }
        }
    }

    let input = InputBuffer::from_bytes(source.into_bytes());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(parsed) => parsed,
        Err(findings) => panic!("hbond fixture failed: {findings:?}"),
    };
    let count = structure.atom_count();
    let mut data = structure.data().clone();
    let oxygens = |atoms: u32| match AnnotationColumn::from_entries((0..atoms).map(|atom| {
        if atom % 2 == 0 {
            (true, Presence::Present)
        } else {
            (false, Presence::Inapplicable)
        }
    })) {
        Ok(column) => column,
        Err(error) => panic!("hbond fixture annotation failed: {error:?}"),
    };
    data.annotations.insert(
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(oxygens(count)),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(oxygens(count)),
    );
    let mut bonds = BondTableBuilder::new();
    for oxygen in (0..count).step_by(2) {
        bonds.push(BondRecord {
            atom_a: molframe_core::AtomIndex::new(oxygen),
            atom_b: molframe_core::AtomIndex::new(oxygen + 1),
            order: BondOrder::Single,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = bonds.finish();
    Structure::new(data)
}

/// Hydrogen-bond detection on a 20³ water lattice (40k atoms). The periodic
/// brute-force row is the evidence base for routing periodic search through
/// the cost model.
pub(super) fn bench_hydrogen_bonds(c: &mut Criterion) {
    let structure = hydrogen_bond_lattice(20);
    let cell = molframe_core::structure::UnitCell {
        lengths: [56.0, 56.0, 56.0],
        angles: [90.0; 3],
    };
    let mut data = structure.data().clone();
    data.cell = Some(cell);
    let periodic = Structure::new(data);
    let mut group = c.benchmark_group("analysis_hbond");
    group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
    for (label, structure, periodic) in [
        ("cell_list", &structure, false),
        ("cell_list/periodic", &periodic, true),
        ("brute_force/periodic", &periodic, true),
    ] {
        let backend = if label.starts_with("brute_force") {
            SpatialBackend::BruteForce
        } else {
            SpatialBackend::CellList
        };
        group.bench_function(label, |b| {
            b.iter(|| {
                let bonds = match hydrogen_bonds(
                    structure,
                    HydrogenBondOptions {
                        maximum_donor_acceptor_distance: 3.5,
                        minimum_angle_degrees: 150.0,
                        backend,
                        periodic,
                    },
                    &context(),
                ) {
                    Ok(bonds) => bonds,
                    Err(error) => panic!("hbond benchmark failed: {error:?}"),
                };
                black_box(bonds.len());
            });
        });
    }
    group.finish();
}
