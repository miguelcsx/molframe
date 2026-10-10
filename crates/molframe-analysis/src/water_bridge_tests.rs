use super::{
    WaterBridgeConfidence, WaterBridgeMode, WaterBridgeOptions, WaterBridgePolicy, water_bridges,
    water_bridges_with_policy,
};
use crate::HydrogenBondOptions;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{
    AnnotationColumn, AtomAnnotation, BondOrder, BondProvenance, BondRecord, BondTableBuilder,
    ExecutionContext, Presence,
};
use molframe_spatial::SpatialBackend;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N DON A 1 0 0 0\nATOM 2 H H DON A 1 1 0 0\n\
HETATM 3 O O HOH W 2 2.8 0 0\nHETATM 4 H H1 HOH W 2 3.8 0 0\n\
ATOM 5 O O ACC B 3 5.6 0 0\n";

#[test]
fn two_oriented_hydrogen_bonds_form_one_water_bridge() {
    let bridges = water_bridges(
        &annotated_structure(),
        WaterBridgeOptions {
            hydrogen_bonds: HydrogenBondOptions {
                maximum_donor_acceptor_distance: 3.5,
                minimum_angle_degrees: 150.0,
                backend: SpatialBackend::BruteForce,
                periodic: false,
            },
        },
        &ExecutionContext::default(),
    )
    .expect("valid bridge network");
    assert_eq!(bridges.len(), 1);
    let bridge = bridges.row(0).expect("one bridge");
    assert_eq!(bridge.water.get(), 2);
    assert_eq!(bridge.first.get(), 0);
    assert_eq!(bridge.second.get(), 4);
}

fn annotated_structure() -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(5, &[0, 2])),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(5, &[2, 4])),
    );
    data.annotations.insert(
        molframe_core::COMPONENT_KIND_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_entries((0..5).map(|atom| {
                let kind = if [2, 3].contains(&atom) {
                    molframe_chem::ComponentKind::Solvent
                } else {
                    molframe_chem::ComponentKind::AminoAcid
                };
                (kind.code(), Presence::Present)
            }))
            .expect("small annotation column"),
        ),
    );
    let mut bonds = BondTableBuilder::new();
    for (first, second) in [(0, 1), (2, 3)] {
        bonds.push(BondRecord {
            atom_a: molframe_core::AtomIndex::new(first),
            atom_b: molframe_core::AtomIndex::new(second),
            order: BondOrder::Single,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = bonds.finish();
    molframe_core::Structure::new(data)
}

fn boolean_roles(atoms: u32, selected: &[u32]) -> AnnotationColumn<bool> {
    AnnotationColumn::from_entries((0..atoms).map(|atom| {
        if selected.contains(&atom) {
            (true, Presence::Present)
        } else {
            (false, Presence::Inapplicable)
        }
    }))
    .expect("small annotation column")
}

const ALT_HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn options() -> WaterBridgeOptions {
    WaterBridgeOptions {
        hydrogen_bonds: HydrogenBondOptions {
            maximum_donor_acceptor_distance: 3.5,
            minimum_angle_degrees: 150.0,
            backend: SpatialBackend::BruteForce,
            periodic: false,
        },
    }
}

/// Builds a structure whose `solvent` atoms carry the solvent kind. Bonds join
/// each listed pair; every atom listed in `donors` or `acceptors` has the role.
fn custom(
    source: &str,
    solvent: &[u32],
    donors: &[u32],
    acceptors: &[u32],
    bonds: &[(u32, u32)],
) -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let count = structure.atom_count();
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(count, donors)),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(count, acceptors)),
    );
    data.annotations.insert(
        molframe_core::COMPONENT_KIND_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_entries((0..count).map(|atom| {
                let kind = if solvent.contains(&atom) {
                    molframe_chem::ComponentKind::Solvent
                } else {
                    molframe_chem::ComponentKind::AminoAcid
                };
                (kind.code(), Presence::Present)
            }))
            .expect("small annotation column"),
        ),
    );
    let mut table = BondTableBuilder::new();
    for (first, second) in bonds {
        table.push(BondRecord {
            atom_a: molframe_core::AtomIndex::new(*first),
            atom_b: molframe_core::AtomIndex::new(*second),
            order: BondOrder::Single,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = table.finish();
    molframe_core::Structure::new(data)
}

fn bridge_count(structure: &molframe_core::Structure) -> usize {
    water_bridges(structure, options(), &ExecutionContext::default())
        .expect("valid")
        .len()
}

/// N-H ... O(solvent)-H ... O across three residues; `solvent_comp` names the
/// middle residue and `extra` adds atoms to it.
fn solvent_source(solvent_comp: &str, extra: &str) -> String {
    format!(
        "{ALT_HEADER}ATOM 1 N N . DON A 1 0 0 0\nATOM 2 H H . DON A 1 1 0 0\n\
HETATM 3 O O . {solvent_comp} W 2 2.8 0 0\nHETATM 4 H H1 . {solvent_comp} W 2 3.8 0 0\n\
{extra}ATOM 99 O O . ACC B 3 5.6 0 0\n"
    )
}

/// `extra` is the number of atoms `solvent_source` inserted after the water.
fn solvent_structure(source: &str, extra: u32) -> molframe_core::Structure {
    let solvent: Vec<u32> = (2..4 + extra).collect();
    custom(
        source,
        &solvent,
        &[0, 2],
        &[2, 4 + extra],
        &[(0, 1), (2, 3)],
    )
}

#[test]
fn a_non_water_solvent_does_not_bridge() {
    let glycerol = "HETATM 5 C C1 . GOL W 2 2.8 1.5 0\n";
    let gol = solvent_structure(&solvent_source("GOL", glycerol), 1);
    assert_eq!(bridge_count(&gol), 0);
    // The same geometry as water, with no carbon, does bridge.
    assert_eq!(
        bridge_count(&solvent_structure(&solvent_source("HOH", ""), 0)),
        1
    );
}

#[test]
fn water_is_also_recognised_by_composition() {
    // An unfamiliar component id with one oxygen and two hydrogens is water.
    let extra = "HETATM 6 H H2 . XYZ W 2 2.8 1 0\n";
    let structure = custom(
        &solvent_source("XYZ", extra),
        &[2, 3, 4],
        &[0, 2],
        &[2, 5],
        &[(0, 1), (2, 3)],
    );
    assert_eq!(bridge_count(&structure), 1);
}

#[test]
fn partners_in_incompatible_alternate_locations_do_not_bridge() {
    let source = |acceptor_alt: &str, donor_alt: &str| {
        format!(
            "{ALT_HEADER}ATOM 1 N N {donor_alt} DON A 1 0 0 0\nATOM 2 H H {donor_alt} DON A 1 1 0 0\n\
HETATM 3 O O . HOH W 2 2.8 0 0\nHETATM 4 H H1 . HOH W 2 3.8 0 0\n\
ATOM 5 O O {acceptor_alt} ACC B 3 5.6 0 0\n"
        )
    };
    assert_eq!(bridge_count(&solvent_structure(&source("B", "A"), 0)), 0);
    assert_eq!(bridge_count(&solvent_structure(&source("A", "A"), 0)), 1);
}

#[test]
fn heavy_atom_mode_bridges_hydrogen_free_waters_and_flags_it() {
    let source = format!(
        "{ALT_HEADER}ATOM 1 N N . DON A 1 0 0 0\nATOM 2 H H . DON A 1 1 0 0\n\
HETATM 3 O O . HOH W 2 2.8 0 0\nATOM 4 O O . ACC B 3 5.6 0 0\n"
    );
    let structure = custom(&source, &[2], &[0, 2], &[2, 3], &[(0, 1)]);
    let context = ExecutionContext::default();

    let strict = water_bridges_with_policy(
        &structure,
        options(),
        WaterBridgePolicy::default(),
        &context,
    )
    .expect("strict");
    assert!(strict.bridges.is_empty());
    assert_eq!(strict.confidence, WaterBridgeConfidence::Hydrogen);
    assert_eq!(strict.waters_without_hydrogens, 1);

    let policy = WaterBridgePolicy {
        mode: WaterBridgeMode::HeavyAtomOnly,
        ..WaterBridgePolicy::default()
    };
    let heavy = water_bridges_with_policy(&structure, options(), policy, &context).expect("heavy");
    assert_eq!(heavy.bridges.len(), 1);
    assert_eq!(heavy.confidence, WaterBridgeConfidence::HeavyAtomOnly);
}
