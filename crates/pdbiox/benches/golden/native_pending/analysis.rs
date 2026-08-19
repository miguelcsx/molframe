use std::fmt::Debug;

use criterion::{BenchmarkGroup, Throughput, black_box};
use pdbiox::{
    AnnotationColumn, AtomAnnotation, AtomIndex, BondOrder, BondProvenance, BondRecord,
    BondTableBuilder, Presence, ReadOptions, Structure,
};

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
    let grid = pdbiox::validate::ReferenceDistribution::grid(
        "general",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![1.0, 1.0, 1.0, 1.0],
    )
    .required("GW-010 grid failed");
    let library = pdbiox::validate::ReferenceLibrary::new("rama", "golden-1", [grid])
        .required("GW-010 library failed");
    let basin = pdbiox::validate::RamachandranBasin::new(
        pdbiox::validate::RamachandranRegion::AlphaHelixRight,
        "general",
    )
    .required("GW-010 basin failed");
    let options = pdbiox::validate::RamachandranOptions::new(&library, [basin], 0.0)
        .required("GW-010 options failed");
    group.throughput(Throughput::Elements(structure.residue_count() as u64));
    group.bench_function("GW-010", |b| {
        b.iter(|| {
            let records = pdbiox::validate::ramachandran(&structure, &options)
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
            let bonds = pdbiox::analysis::hydrogen_bonds(
                &structure,
                pdbiox::analysis::HydrogenBondOptions {
                    maximum_donor_acceptor_distance: 3.5,
                    minimum_angle_degrees: 150.0,
                    backend: pdbiox::SpatialBackend::BruteForce,
                    periodic: false,
                },
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
            let records = pdbiox::analysis::secondary_structure(
                &structure,
                &pdbiox::analysis::DsspOptions {
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
        pdbiox::core::contract::ContactDefinition::DistanceCutoff { tolerance: 0.0 },
        pdbiox::core::contract::ContactDefinition::DistanceCutoff { tolerance: 0.5 },
        pdbiox::core::contract::ContactDefinition::DistanceCutoff { tolerance: 1.0 },
        pdbiox::core::contract::ContactDefinition::SurfaceBased { probe: 1.2 },
        pdbiox::core::contract::ContactDefinition::SurfaceBased { probe: 1.4 },
    ];
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-038", |b| {
        b.iter(|| {
            let fingerprints: Vec<_> = policies
                .into_iter()
                .map(|contact_def| {
                    let policy = pdbiox::AnalysisPolicy {
                        contact_def,
                        ..pdbiox::AnalysisPolicy::default()
                    };
                    let count = match contact_def {
                        pdbiox::core::contract::ContactDefinition::DistanceCutoff { tolerance } => {
                            pdbiox::analysis::atom_contacts(
                                &structure,
                                2.5 + tolerance,
                                pdbiox::SpatialBackend::BruteForce,
                            )
                            .required("GW-038 distance failed")
                            .len()
                        }
                        pdbiox::core::contract::ContactDefinition::SurfaceBased { probe } => {
                            pdbiox::analysis::surface_contacts(
                                &structure,
                                &[1.7, 1.7],
                                0.5,
                                probe,
                                2.0,
                                0.1,
                                pdbiox::SpatialBackend::BruteForce,
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
    let mut data = structure.data().clone();
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
        pdbiox::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(0)),
    );
    data.annotations.insert(
        pdbiox::HBOND_ACCEPTOR_ANNOTATION,
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
    Structure::new(data)
}

fn with_polymer_roles(structure: &Structure) -> Structure {
    let values: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("N") => pdbiox::PolymerAtomRole::PROTEIN_NITROGEN,
                Some("CA") => pdbiox::PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("C") => pdbiox::PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                Some("O") => pdbiox::PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                _ => pdbiox::PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.data().clone();
    data.annotations.insert(
        pdbiox::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_entries(values).required("GW-016 role annotation failed"),
        ),
    );
    Structure::new(data)
}

fn rama_structure() -> Structure {
    let source = read(RAMA_CIF);
    let mut data = source.data().clone();
    let roles = [
        pdbiox::PolymerAtomRole::PROTEIN_NITROGEN,
        pdbiox::PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        pdbiox::PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    let values = AnnotationColumn::from_entries(
        (0..9).map(|index| (roles[index % roles.len()].code(), Presence::Present)),
    )
    .required("GW-010 role annotation failed");
    data.annotations.insert(
        pdbiox::POLYMER_ATOM_ROLE_ANNOTATION,
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
    Structure::new(data)
}

fn read(source: &str) -> Structure {
    pdbiox::read_bytes(
        source.as_bytes().to_vec(),
        Some("golden.cif"),
        &ReadOptions::new(),
    )
    .required("golden benchmark fixture failed")
    .0
}
