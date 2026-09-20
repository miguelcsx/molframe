use std::fmt::Debug;

use criterion::{BenchmarkGroup, Throughput, black_box};
use molframe::{
    AnnotationColumn, AtomAnnotation, AtomIndex, BondOrder, BondProvenance, BondRecord,
    BondTableBuilder, ReadOptions, Structure,
};
use molframe_core::column::Presence;
use molframe_core::structure::Structure as CoreStructure;

trait BenchRequired<T> {
    fn required(self, context: &str) -> T;
}

impl<T, E: Debug> BenchRequired<T> for Result<T, E> {
    fn required(self, context: &str) -> T {
        match self {
            Ok(value) => value,
            Err(error) => panic!("{context}: {error:?}"),
        }
    }
}

impl<T> BenchRequired<T> for Option<T> {
    fn required(self, context: &str) -> T {
        match self {
            Some(value) => value,
            None => panic!("{context}"),
        }
    }
}

const HBOND_CIF: &str = r"data_hbond
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N DON A 1 0 0 0
ATOM 2 H H DON A 1 1 0 0
ATOM 3 O O ACC A 2 2.8 0 0
";

const DSSP_CIF: &str = r"data_dssp
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 0 0 0
ATOM 2 C CA GLY A 1 1.5 0 0
ATOM 3 C C GLY A 1 2 1 0
ATOM 4 O O GLY A 1 3 1 0
ATOM 5 N N GLY A 2 1.5 2 0
ATOM 6 C CA GLY A 2 2 3 0
ATOM 7 C C GLY A 2 3 3 0
ATOM 8 O O GLY A 2 4 3 0
";

const RAMA_CIF: &str = r"data_rama
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N ALA A 1 0 0 0
ATOM 2 C CA ALA A 1 1.46 0 0
ATOM 3 C C ALA A 1 2.0 1.2 0
ATOM 4 N N ALA A 2 3.3 1.4 0
ATOM 5 C CA ALA A 2 4.0 2.5 0
ATOM 6 C C ALA A 2 5.4 2.5 0
ATOM 7 N N ALA A 3 6.0 3.6 0
ATOM 8 C CA ALA A 3 7.0 3.6 0
ATOM 9 C C ALA A 3 8.0 4.0 0
";

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    bench_gw_015(group);
    bench_gw_016(group);
    bench_gw_038(group);
    bench_gw_010(group);
}

fn bench_gw_010(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = rama_structure();
    let grid = molframe::validation::ReferenceDistribution::grid(
        "general",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![1.0, 1.0, 1.0, 1.0],
    )
    .required("GW-010 grid failed");
    let library = molframe::validation::ReferenceLibrary::new("rama", "golden-1", [grid])
        .required("GW-010 library failed");
    let basin = molframe::validation::RamachandranBasin::new(
        molframe::validation::RamachandranRegion::AlphaHelixRight,
        "general",
    )
    .required("GW-010 basin failed");
    let options = molframe::validation::RamachandranOptions::new(&library, [basin], 0.0)
        .required("GW-010 options failed");
    group.throughput(Throughput::Elements(structure.residue_count() as u64));
    group.bench_function("GW-010", |b| {
        b.iter(|| {
            let records = molframe::validation::ramachandran(structure.engine(), &options)
                .required("GW-010 classification failed");
            black_box(records.len());
        });
    });
}

fn bench_gw_015(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = annotated_hbond_structure();
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-015", |b| {
        b.iter(|| {
            let bonds = molframe::analysis::hydrogen_bonds(
                structure.engine(),
                molframe::analysis::HydrogenBondOptions {
                    maximum_donor_acceptor_distance: 3.5,
                    minimum_angle_degrees: 150.0,
                    backend: molframe::spatial::SpatialBackend::BruteForce,
                    periodic: false,
                },
                &molframe::ExecutionContext::default(),
            )
            .required("GW-015 failed");
            black_box(bonds.len());
        });
    });
}

fn bench_gw_016(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let source = read(DSSP_CIF);
    let structure = with_polymer_roles(&source);
    group.throughput(Throughput::Elements(structure.residue_count() as u64));
    group.bench_function("GW-016", |b| {
        b.iter(|| {
            let records = molframe::analysis::secondary_structure(
                structure.engine(),
                &molframe::analysis::DsspOptions {
                    electrostatic_prefactor: 332.0 * 0.42 * 0.20,
                    hydrogen_bond_energy: -0.5,
                    amide_hydrogen_distance: 1.0,
                    minimum_sequence_separation: 2,
                    helix_offset: 4,
                    turn_offsets: 3..=5,
                },
            )
            .required("GW-016 failed");
            black_box(records.len());
        });
    });
}

fn bench_gw_038(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = read(
        r"data_interface
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA ALA A 1 0 0 0
ATOM 2 C CA ALA B 1 2.5 0 0
",
    );
    let policies = [
        molframe_core::contract::ContactDefinition::DistanceCutoff { tolerance: 0.0 },
        molframe_core::contract::ContactDefinition::DistanceCutoff { tolerance: 0.5 },
        molframe_core::contract::ContactDefinition::DistanceCutoff { tolerance: 1.0 },
        molframe_core::contract::ContactDefinition::SurfaceBased { probe: 1.2 },
        molframe_core::contract::ContactDefinition::SurfaceBased { probe: 1.4 },
    ];
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-038", |b| {
        b.iter(|| {
            let fingerprints: Vec<_> = policies
                .into_iter()
                .map(|contact_def| {
                    let policy = molframe::AnalysisPolicy {
                        contact_def,
                        ..molframe::AnalysisPolicy::default()
                    };
                    let count = match contact_def {
                        molframe_core::contract::ContactDefinition::DistanceCutoff {
                            tolerance,
                        } => molframe::analysis::atom_contacts(
                            structure.engine(),
                            2.5 + tolerance,
                            molframe::spatial::SpatialBackend::BruteForce,
                            &molframe::ExecutionContext::default(),
                        )
                        .required("GW-038 distance failed")
                        .len(),
                        molframe_core::contract::ContactDefinition::SurfaceBased { probe } => {
                            molframe::analysis::surface_contacts(
                                structure.engine(),
                                &[1.7, 1.7],
                                molframe::analysis::SurfaceContactOptions {
                                    tolerance: 0.5,
                                    probe,
                                    surface_density: 2.0,
                                    minimum_area: 0.1,
                                    backend: molframe::spatial::SpatialBackend::BruteForce,
                                },
                                &molframe::ExecutionContext::default(),
                            )
                            .required("GW-038 surface failed")
                            .len()
                        }
                        _ => 0,
                    };
                    (policy.fingerprint(), count)
                })
                .collect();
            black_box(fingerprints);
        });
    });
}

fn annotated_hbond_structure() -> Structure {
    let structure = read(HBOND_CIF);
    let mut data = structure.engine().data().clone();
    let roles = |selected: u32| {
        AnnotationColumn::from_entries((0..3).map(|atom| {
            if atom == selected {
                (true, Presence::Present)
            } else {
                (false, Presence::Inapplicable)
            }
        }))
        .required("GW-015 annotation failed")
    };
    data.annotations.insert(
        molframe::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(0)),
    );
    data.annotations.insert(
        molframe::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(2)),
    );
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(0),
        atom_b: AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::ChemicalComponentDictionary,
    });
    data.bonds = bonds.finish();
    CoreStructure::new(data).into()
}

fn with_polymer_roles(structure: &Structure) -> Structure {
    let values: Vec<_> = structure
        .engine()
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("N") => molframe::chemistry::PolymerAtomRole::PROTEIN_NITROGEN,
                Some("CA") => molframe::chemistry::PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("C") => molframe::chemistry::PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                Some("O") => molframe::chemistry::PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                _ => molframe::chemistry::PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.engine().data().clone();
    data.annotations.insert(
        molframe::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_entries(values).required("GW-016 role annotation failed"),
        ),
    );
    CoreStructure::new(data).into()
}

fn rama_structure() -> Structure {
    let source = read(RAMA_CIF);
    let mut data = source.engine().data().clone();
    let roles = [
        molframe::chemistry::PolymerAtomRole::PROTEIN_NITROGEN,
        molframe::chemistry::PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        molframe::chemistry::PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    let values = AnnotationColumn::from_entries(
        (0..9).map(|index| (roles[index % roles.len()].code(), Presence::Present)),
    )
    .required("GW-010 role annotation failed");
    data.annotations.insert(
        molframe::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(values),
    );
    let mut bonds = BondTableBuilder::new();
    for (carbon, nitrogen) in [(2, 3), (5, 6)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(carbon),
            atom_b: AtomIndex::new(nitrogen),
            order: BondOrder::Single,
            provenance: BondProvenance::User,
        });
    }
    data.bonds = bonds.finish();
    CoreStructure::new(data).into()
}

fn read(source: &str) -> Structure {
    molframe::read_bytes(
        source.as_bytes().to_vec(),
        Some("golden.cif"),
        &ReadOptions::new(),
    )
    .required("golden benchmark fixture failed")
    .0
}
