use super::{
    HydrogenBond, HydrogenBondError, HydrogenBondOptions, hydrogen_bonds as hydrogen_bonds_native,
};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{
    AnnotationColumn, AtomAnnotation, BondOrder, BondProvenance, BondRecord, BondTableBuilder,
    ExecutionContext, Presence,
};
use molframe_spatial::SpatialBackend;

fn hydrogen_bonds(
    structure: &molframe_core::Structure,
    options: HydrogenBondOptions,
) -> Result<Vec<HydrogenBond>, HydrogenBondError> {
    hydrogen_bonds_native(structure, options, &ExecutionContext::default())
}

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
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
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
        backend: SpatialBackend::CellList,
        periodic: false,
    }
}

fn annotated_structure(source: &str) -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
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
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(0)),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(roles(2)),
    );
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: molframe_core::AtomIndex::new(0),
        atom_b: molframe_core::AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::ChemicalComponentDictionary,
    });
    data.bonds = bonds.finish();
    molframe_core::Structure::new(data)
}

/// A lattice of water-like O-H pairs, every oxygen both donor and acceptor.
///
/// Large enough to span several query blocks, so a parallel search genuinely
/// divides the work rather than collapsing to one accumulator.
fn water_lattice(side: i16) -> molframe_core::Structure {
    use std::fmt::Write;

    let mut source = String::from(
        "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut serial_id = 0u32;
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                let (ox, oy, oz) = (f32::from(x) * 2.8, f32::from(y) * 2.8, f32::from(z) * 2.8);
                let residue = serial_id / 2 + 1;
                serial_id += 1;
                writeln!(
                    source,
                    "ATOM {serial_id} O O HOH A {residue} {ox:.3} {oy:.3} {oz:.3}"
                )
                .expect("fixture row");
                serial_id += 1;
                writeln!(
                    source,
                    "ATOM {serial_id} H H1 HOH A {residue} {:.3} {oy:.3} {oz:.3}",
                    ox + 0.96
                )
                .expect("fixture row");
            }
        }
    }

    let input = InputBuffer::from_bytes(source.into_bytes());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let count = structure.atom_count();
    let mut data = structure.data().clone();

    // Even indices are the oxygens, which both donate and accept.
    let oxygens = || {
        AnnotationColumn::from_entries((0..count).map(|atom| {
            if atom % 2 == 0 {
                (true, Presence::Present)
            } else {
                (false, Presence::Inapplicable)
            }
        }))
        .expect("annotation column")
    };
    data.annotations.insert(
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(oxygens()),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(oxygens()),
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
    molframe_core::Structure::new(data)
}

#[test]
fn the_bond_list_is_identical_at_every_worker_count() {
    let structure = water_lattice(5);
    let settings = options(90.0);
    let serial = hydrogen_bonds(&structure, settings).expect("serial search");
    assert!(!serial.is_empty(), "the fixture must produce bonds");

    for workers in [2, 4, 8, 16] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let parallel =
            hydrogen_bonds_native(&structure, settings, &context).expect("parallel search");
        assert_eq!(
            parallel, serial,
            "worker count {workers} changed the result"
        );
    }
}
