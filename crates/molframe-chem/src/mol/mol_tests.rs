use super::{MolVersion, parse_mol_record, parse_sdf_records, write_mol, write_sdf};
use molframe_core::element::Element;

const CARBON_MONOXIDE: &str = "CO
  molframe
comment
  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.2300    0.0000    0.0000 O   0  0
  1  2  2  0
M  CHG  1   2  -1
M  END
>  <SOURCE>
fixture

";

#[test]
fn v2000_roundtrip_preserves_headers_charge_and_property() {
    let record = parse_mol_record(CARBON_MONOXIDE).expect("valid V2000");
    assert_eq!(record.atom_metadata[1].formal_charge, Some(-1));
    let encoded = write_mol(&record).expect("write V2000");
    let decoded = parse_mol_record(&encoded).expect("read written V2000");
    assert_eq!(decoded, record);
}

#[test]
fn v3000_roundtrip_preserves_atom_and_bond_attributes() {
    let source = "CO\nmolframe\ncomment\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n\
M  V30 7 C 0 0 0 0 MASS=13\nM  V30 9 O 1.2 0 0 0 CHG=-1\nM  V30 END ATOM\n\
M  V30 BEGIN BOND\nM  V30 1 2 7 9 CFG=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    let record = parse_mol_record(source).expect("valid V3000");
    assert_eq!(record.version, MolVersion::V3000);
    let encoded = write_mol(&record).expect("write V3000");
    let decoded = parse_mol_record(&encoded).expect("read written V3000");
    assert_eq!(decoded, record);
}

#[test]
fn sdf_writes_and_reads_multiple_records_in_order() {
    let first = parse_mol_record(CARBON_MONOXIDE).expect("first record");
    let mut second = first.clone();
    second.name = "second".into();
    let encoded = write_sdf(&[first.clone(), second.clone()]).expect("write SDF");
    assert_eq!(
        parse_sdf_records(&encoded).expect("read SDF"),
        vec![first, second]
    );
}

#[test]
fn complete_parser_exposes_the_graph_without_discarding_metadata() {
    let record = parse_mol_record(CARBON_MONOXIDE).expect("valid MOL block");
    assert_eq!(record.molecule.atoms[0].element, Element::CARBON);
    assert_eq!(record.molecule.bonds[0].order, 2);
}
