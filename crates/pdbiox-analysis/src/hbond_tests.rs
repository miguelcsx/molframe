use super::{HydrogenBondError, HydrogenBondOptions, hydrogen_bonds};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{
    AnnotationColumn, AtomAnnotation, BondOrder, BondProvenance, BondRecord, BondTableBuilder,
    Presence,
};
use pdbiox_spatial::SpatialBackend;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N DON A 1 0 0 0\n\
ATOM 2 H H DON A 1 1 0 0\n\
ATOM 3 O O ACC A 2 2.8 0 0\n";

#[test]
fn explicit_hydrogen_and_ccd_roles_define_direction_and_angle() {
    let structure = annotated_structure(SOURCE);
    let bonds = hydrogen_bonds(&structure, options(150.0)).expect("valid hydrogen bonds");
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].donor.get(), 0);
    assert_eq!(bonds[0].hydrogen.get(), 1);
    assert_eq!(bonds[0].acceptor.get(), 2);
    assert!((bonds[0].donor_acceptor_distance - 2.8).abs() < 1.0e-5);
    assert!((bonds[0].angle_degrees - 180.0).abs() < 1.0e-4);
}

#[test]
fn angular_policy_filters_a_bent_geometry() {
    let bent = SOURCE.replace("ATOM 3 O O ACC A 2 2.8 0 0", "ATOM 3 O O ACC A 2 1 1.8 0");
    let structure = annotated_structure(&bent);
    assert!(
        hydrogen_bonds(&structure, options(150.0))
            .expect("valid request")
            .is_empty()
    );
}

#[test]
fn absent_ccd_roles_are_an_error_not_a_heavy_atom_fallback() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    assert!(matches!(
        hydrogen_bonds(&structure, options(150.0)),
        Err(HydrogenBondError::MissingChemistry)
    ));
}

fn options(minimum_angle_degrees: f64) -> HydrogenBondOptions {
    HydrogenBondOptions {
        maximum_donor_acceptor_distance: 3.5,
        minimum_angle_degrees,
        backend: SpatialBackend::BruteForce,
        periodic: false,
    }
}

fn annotated_structure(source: &str) -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let mut data = structure.data().clone();
    let roles = |selected: u32| {
        AnnotationColumn::from_entries((0..3).map(|atom| {
            if atom == selected {
                (true, Presence::Present)
            } else {
                (false, Presence::Inapplicable)
            }
        }))
        .expect("small annotation column")
    };
    data.annotations.insert(
        pdbiox_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(0)),
    );
    data.annotations.insert(
        pdbiox_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(2)),
    );
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: pdbiox_core::AtomIndex::new(0),
        atom_b: pdbiox_core::AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::ChemicalComponentDictionary,
    });
    data.bonds = bonds.finish();
    pdbiox_core::Structure::new(data)
}
