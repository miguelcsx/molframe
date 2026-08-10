use super::*;
use pdbiox_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AtomAnnotation, PARTIAL_CHARGE_ANNOTATION,
};
use pdbiox_core::{Format, InputBuffer, ReadOptions};

const PQR: &str = "ATOM      1  N   ALA A   1      11.104  13.207   9.124 -0.3000 1.5500\nEND\n";
const PDBQT: &str = "ROOT\n\
ATOM      1  N   LIG A   1      11.104  13.207   9.124  1.00  0.00      -0.300 N \n\
ENDROOT\nTORSDOF 0\n";

#[test]
fn pqr_charge_and_radius_survive_a_canonical_round_trip() {
    let structure = read_variant(PQR, Format::Pqr);
    assert_real(&structure, PARTIAL_CHARGE_ANNOTATION, -0.3);
    assert_real(&structure, ATOM_RADIUS_ANNOTATION, 1.55);
    let written = match write_pqr(&structure, &PdbOptions::new()) {
        Ok(written) => written,
        Err(findings) => panic!("PQR write failed: {findings:?}"),
    };
    let round_tripped = read_variant(&written, Format::Pqr);
    assert_real(&round_tripped, PARTIAL_CHARGE_ANNOTATION, -0.3);
    assert_real(&round_tripped, ATOM_RADIUS_ANNOTATION, 1.55);
}

#[test]
fn pdbqt_charge_and_type_survive_a_rigid_round_trip() {
    let structure = read_variant(PDBQT, Format::Pdbqt);
    assert_real(&structure, PARTIAL_CHARGE_ANNOTATION, -0.3);
    let Some(AtomAnnotation::Symbol(types)) = structure.annotations().get(AUTODOCK_TYPE_ANNOTATION)
    else {
        panic!("AutoDock types absent")
    };
    let Some((atom_type, _)) = types.get(0) else {
        panic!("first AutoDock type absent")
    };
    assert_eq!(structure.resolve(atom_type), Some("N"));
    let written = match write_pdbqt(&structure, &PdbOptions::new()) {
        Ok(written) => written,
        Err(findings) => panic!("PDBQT write failed: {findings:?}"),
    };
    assert!(written.starts_with("ROOT\n"));
    assert!(written.ends_with("TORSDOF 0\n"));
    let round_tripped = read_variant(&written, Format::Pdbqt);
    assert_real(&round_tripped, PARTIAL_CHARGE_ANNOTATION, -0.3);
}

fn read_variant(text: &str, format: Format) -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().format(format);
    let result = match format {
        Format::Pqr => crate::read_pqr(&input, &options),
        Format::Pdbqt => crate::read_pdbqt(&input, &options),
        _ => panic!("test format is not a charged PDB variant"),
    };
    match result {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("variant read failed: {findings:?}"),
    }
}

fn assert_real(structure: &pdbiox_core::Structure, name: &str, expected: f64) {
    let Some(AtomAnnotation::Real(column)) = structure.annotations().get(name) else {
        panic!("{name} absent")
    };
    let Some((value, _)) = column.get(0) else {
        panic!("{name} row absent")
    };
    assert!((value - expected).abs() < 1.0e-8);
}
