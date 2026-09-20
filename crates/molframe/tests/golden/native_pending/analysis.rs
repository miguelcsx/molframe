use molframe::{
    AnnotationColumn, AtomAnnotation, AtomIndex, BondOrder, BondProvenance, BondRecord,
    BondTableBuilder, ReadOptions, Structure,
};
use molframe_core::column::Presence;
use molframe_core::structure::Structure as CoreStructure;

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

#[test]
fn gw_015_detects_oriented_hydrogen_bonds_from_explicit_chemistry() {
    let structure = annotated_hbond_structure();
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
    .unwrap_or_else(|error| panic!("hydrogen-bond workflow failed: {error}"));
    assert_eq!(bonds.len(), 1);
    let bond = bonds.row(0).unwrap_or_else(|| panic!("one hydrogen bond"));
    assert_eq!(
        (bond.donor.get(), bond.hydrogen.get(), bond.acceptor.get()),
        (0, 1, 2)
    );
    assert!((bond.angle_degrees - 180.0).abs() < 1.0e-4);
}

#[test]
fn gw_016_assigns_one_secondary_structure_record_per_backbone_residue() {
    let source = read(DSSP_CIF);
    let structure = with_polymer_roles(&source);
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
    .unwrap_or_else(|error| panic!("secondary-structure workflow failed: {error}"));
    assert_eq!(records.len(), 2);
    assert!(
        records
            .iter()
            .all(|record| record.kind == molframe::analysis::SseKind::Coil)
    );
}

#[test]
fn gw_038_reports_policy_sensitive_interface_contact_runs() {
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
    let mut fingerprints = Vec::new();
    let mut counts = Vec::new();
    for contact_def in policies {
        let policy = molframe::AnalysisPolicy {
            contact_def,
            ..molframe::AnalysisPolicy::default()
        };
        let count = match contact_def {
            molframe_core::contract::ContactDefinition::DistanceCutoff { tolerance } => {
                molframe::analysis::atom_contacts(
                    structure.engine(),
                    2.5 + tolerance,
                    molframe::spatial::SpatialBackend::BruteForce,
                    &molframe::ExecutionContext::default(),
                )
                .unwrap_or_else(|error| panic!("distance contact workflow failed: {error}"))
                .len()
            }
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
                .unwrap_or_else(|error| panic!("surface contact workflow failed: {error}"))
                .len()
            }
            _ => 0,
        };
        fingerprints.push(policy.fingerprint());
        counts.push(count);
    }
    assert_eq!(fingerprints.len(), 5);
    assert!(fingerprints.windows(2).all(|pair| pair[0] != pair[1]));
    assert!(counts.iter().any(|count| *count > 0));
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
        .unwrap_or_else(|error| panic!("annotation fixture failed: {error}"))
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
            AnnotationColumn::from_entries(values)
                .unwrap_or_else(|error| panic!("role fixture failed: {error}")),
        ),
    );
    CoreStructure::new(data).into()
}

fn read(source: &str) -> Structure {
    match molframe::read_bytes(
        source.as_bytes().to_vec(),
        Some("golden.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("golden native fixture failed: {findings:?}"),
    }
}
