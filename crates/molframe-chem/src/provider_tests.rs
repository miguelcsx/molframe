use super::*;
use molframe_core::io::InputBuffer;

const COMPONENT: &str = "data_GLY\n\
_chem_comp.id GLY\n_chem_comp.name GLYCINE\n\
_chem_comp.type 'L-peptide linking'\n_chem_comp.one_letter_code G\n_chem_comp.formula 'C2 H5 N O2'\n\
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
    let (document, _) = molframe_cif::parse(&input).expect("CCD parses");
    let provider = CifProvider::from_document(&document, DictionaryVersion::new("2026-08-01"))
        .expect("component lowers");
    let component = provider
        .get("GLY")
        .expect("provider works")
        .expect("GLY exists");
    assert_eq!(component.kind, ComponentKind::AminoAcid);
    assert_eq!(component.one_letter_code, Some(b'G'));
    assert_eq!(component.atoms.len(), 2);
    assert_eq!(component.bonds[0].order, BondOrder::Single);
    assert_eq!(provider.version().as_str(), "2026-08-01");
}

#[test]
fn an_unknown_component_is_absence_not_a_fabricated_definition() {
    let provider =
        MemoryProvider::new(DictionaryVersion::new("test"), []).expect("empty fixture is unique");
    assert!(provider.get("XXX").expect("lookup works").is_none());
}

#[test]
fn in_memory_provider_rejects_duplicate_component_identifiers() {
    let component = Component {
        id: "DUP".into(),
        name: "duplicate".into(),
        kind: ComponentKind::Unknown,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: Arc::from([]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    };
    let result = MemoryProvider::new(
        DictionaryVersion::new("test"),
        [component.clone(), component],
    );
    assert!(matches!(result, Err(ref finding) if finding.code() == Code::E2005));
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

#[test]
fn component_kind_uses_ccd_type_tokens_not_substring_guesses() {
    assert_eq!(kind(Some("L-peptide linking")), ComponentKind::AminoAcid);
    assert_eq!(
        kind(Some("DNA OH 3 prime terminus")),
        ComponentKind::Nucleotide
    );
    assert_eq!(kind(Some("D-saccharide")), ComponentKind::Saccharide);
    assert_eq!(kind(Some("water")), ComponentKind::Solvent);
    assert_eq!(kind(Some("non-polymer")), ComponentKind::NonPolymer);
    assert_eq!(kind(Some("notwater")), ComponentKind::Unknown);
    assert_eq!(kind(Some("polymerase inhibitor")), ComponentKind::Unknown);
}

#[test]
fn monoatomic_charged_non_polymers_are_classified_as_ions_from_chemistry() {
    let atom = ComponentAtom {
        name: "ZN".into(),
        alternate_name: None,
        element: Element::ZINC,
        charge: 2,
        aromatic: false,
        leaving: false,
        stereo: None,
    };
    assert_eq!(
        classify_component("NON-POLYMER", &[atom], &[]),
        ComponentKind::Ion
    );
}

#[test]
fn missing_atom_chemistry_is_rejected_instead_of_becoming_neutral() {
    let source = COMPONENT.replace("N N 0 N N N", "N N ? N N N");
    let input = InputBuffer::from_bytes(source.into_bytes());
    let (document, _) = molframe_cif::parse(&input).expect("CCD syntax parses");
    let findings = CifProvider::from_document(&document, DictionaryVersion::new("test"))
        .expect_err("unknown charge must not become zero");
    assert!(findings.iter().any(|finding| finding.code() == Code::E2002));
}

#[test]
fn invalid_chemical_flags_are_not_interpreted_as_false() {
    let source = COMPONENT.replace("N N 0 N N N", "N N 0 MAYBE N N");
    let input = InputBuffer::from_bytes(source.into_bytes());
    let (document, _) = molframe_cif::parse(&input).expect("CCD syntax parses");
    let findings = CifProvider::from_document(&document, DictionaryVersion::new("test"))
        .expect_err("invalid aromaticity must be reported");
    assert!(findings.iter().any(|finding| finding.code() == Code::E2003));
}

#[test]
fn duplicate_component_blocks_are_rejected_instead_of_replacing_data() {
    let input = InputBuffer::from_bytes(format!("{COMPONENT}\n{COMPONENT}").into_bytes());
    let (document, _) = molframe_cif::parse(&input).expect("CCD syntax parses");
    let findings = CifProvider::from_document(&document, DictionaryVersion::new("test"))
        .expect_err("duplicate identifiers must be rejected");
    assert!(findings.iter().any(|finding| finding.code() == Code::E2005));
}

#[test]
fn unresolved_ccd_bond_endpoints_are_rejected() {
    let source = COMPONENT.replace("N CA SING N N", "N ZZ SING N N");
    let input = InputBuffer::from_bytes(source.into_bytes());
    let (document, _) = molframe_cif::parse(&input).expect("CCD syntax parses");
    let findings = CifProvider::from_document(&document, DictionaryVersion::new("test"))
        .expect_err("dangling chemical bonds must be rejected");
    assert!(findings.iter().any(|finding| finding.code() == Code::E2005));
}
