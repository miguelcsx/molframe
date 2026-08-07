use super::*;
use pdbiox_core::io::InputBuffer;

const COMPONENT: &str = "data_GLY\n\
_chem_comp.id GLY\n_chem_comp.name GLYCINE\n\
_chem_comp.type 'L-peptide linking'\n_chem_comp.formula 'C2 H5 N O2'\n\
loop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
_chem_comp_atom.pdbx_leaving_atom_flag\n_chem_comp_atom.pdbx_stereo_config\n\
N N 0 N N N\nCA C 0 N N S\n\
loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n\
_chem_comp_bond.pdbx_stereo_config\nN CA SING N N\n";

#[test]
fn a_ccd_document_becomes_a_versioned_component_provider() {
    let input = InputBuffer::from_bytes(COMPONENT.as_bytes().to_vec());
    let (document, _) = pdbiox_cif::parse(&input).expect("CCD parses");
    let provider = CifProvider::from_document(&document, DictionaryVersion::new("2026-08-01"))
        .expect("component lowers");
    let component = provider
        .get("GLY")
        .expect("provider works")
        .expect("GLY exists");
    assert_eq!(component.kind, ComponentKind::AminoAcid);
    assert_eq!(component.atoms.len(), 2);
    assert_eq!(component.bonds[0].order, BondOrder::Single);
    assert_eq!(provider.version().as_str(), "2026-08-01");
}

#[test]
fn an_unknown_component_is_absence_not_a_fabricated_definition() {
    let provider = MemoryProvider::new(DictionaryVersion::new("test"), []);
    assert!(provider.get("XXX").expect("lookup works").is_none());
}

#[test]
fn ccd_input_uses_the_lossless_cif_parser_and_retains_its_version() {
    let input = InputBuffer::from_bytes(COMPONENT.as_bytes().to_vec());
    let (provider, findings) = match read_ccd(&input, DictionaryVersion::new("ccd-test")) {
        Ok(result) => result,
        Err(findings) => panic!("CCD read failed: {findings:?}"),
    };
    assert!(findings.is_empty());
    assert_eq!(provider.version().as_str(), "ccd-test");
    assert!(
        provider
            .get("GLY")
            .is_ok_and(|component| component.is_some())
    );
}
