use super::*;

const CIF: &str = "data_comp_LIG\n_chem_comp.id LIG\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\nLIG C1 C\nLIG C2 C\nLIG O1 O\nLIG N1 N\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.type\n_chem_comp_bond.value_dist\n_chem_comp_bond.value_dist_esd\nLIG C1 C2 coval 1.50 0.02\nloop_\n_chem_comp_angle.comp_id\n_chem_comp_angle.atom_id_1\n_chem_comp_angle.atom_id_2\n_chem_comp_angle.atom_id_3\n_chem_comp_angle.value_angle\n_chem_comp_angle.value_angle_esd\nLIG O1 C1 C2 120 2\nloop_\n_chem_comp_tor.comp_id\n_chem_comp_tor.id\n_chem_comp_tor.atom_id_1\n_chem_comp_tor.atom_id_2\n_chem_comp_tor.atom_id_3\n_chem_comp_tor.atom_id_4\n_chem_comp_tor.value_angle\n_chem_comp_tor.value_angle_esd\n_chem_comp_tor.period\nLIG chi O1 C1 C2 N1 180 10 2\nloop_\n_chem_comp_plane_atom.comp_id\n_chem_comp_plane_atom.plane_id\n_chem_comp_plane_atom.atom_id\n_chem_comp_plane_atom.dist_esd\nLIG p C1 0.02\nLIG p C2 0.02\nLIG p O1 0.02\nloop_\n_chem_comp_chir.comp_id\n_chem_comp_chir.id\n_chem_comp_chir.atom_id_centre\n_chem_comp_chir.atom_id_1\n_chem_comp_chir.atom_id_2\n_chem_comp_chir.atom_id_3\n_chem_comp_chir.volume_sign\nLIG chir C1 C2 O1 N1 positive\n";

fn library(text: &str) -> Result<MonomerLibrary, RestraintError> {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, findings) = molframe_cif::parse(&input).expect("CIF should parse");
    assert!(findings.is_empty());
    lower_monomer_library(&document)
}

#[test]
fn every_restraint_family_is_grouped_and_typed() {
    let library = library(CIF).expect("library should lower");
    let ligand = library.get("LIG").expect("component should exist");
    assert_eq!(ligand.bonds.len(), 1);
    assert_eq!(ligand.angles[0].atoms[1].as_ref(), "C1");
    assert_eq!(ligand.torsions[0].period, 2);
    assert_eq!(ligand.planes[0].atoms.len(), 3);
    assert_eq!(ligand.chirals[0].sign, ChiralVolumeSign::Positive);
}

#[test]
fn absent_atom_references_are_refused() {
    assert_eq!(
        library(&CIF.replace("LIG C1 C2 coval", "LIG C1 XX coval")),
        Err(RestraintError::AtomReference)
    );
}

#[test]
fn official_ccp4_component_when_configured() {
    let Some(path) = std::env::var_os("MOLFRAME_MONOMER_FIXTURE") else {
        return;
    };
    let input = InputBuffer::from_bytes(std::fs::read(path).expect("fixture should be readable"));
    let (library, _) = read_monomer_library(&input).expect("CCP4 monomer should parse");
    let component = library.iter().next().expect("component should exist").1;
    assert!(!component.bonds.is_empty());
    assert!(!component.angles.is_empty());
}
