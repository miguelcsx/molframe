use super::*;

const XML: &str = "<?xml version=\"1.0\"?><PDBx:datablock datablockName=\"X\" xmlns:PDBx=\"http://pdbml.pdb.org/schema/pdbx-v50.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"><PDBx:entryCategory><PDBx:entry id=\"X\"><PDBx:details>A &amp; B</PDBx:details></PDBx:entry></PDBx:entryCategory><PDBx:testCategory><PDBx:test id=\"1\"><PDBx:value>2.5</PDBx:value><PDBx:optional xsi:nil=\"true\"/></PDBx:test><PDBx:test id=\"2\"><PDBx:late>yes</PDBx:late></PDBx:test></PDBx:testCategory></PDBx:datablock>";

#[test]
fn attributes_elements_nil_and_late_columns_become_one_document() {
    let document = parse_pdbml_document(XML.as_bytes()).expect("PDBML should parse");
    let block = document.first_block().expect("block should exist");
    assert_eq!(block.name(), "X");
    let entry = block.category("entry").expect("entry should exist");
    assert_eq!(entry.text("details", 0), Some("A & B"));
    let test = block.category("test").expect("test should exist");
    assert_eq!(
        test.value("value", 0).and_then(CifValue::as_float),
        Some(2.5)
    );
    assert!(matches!(test.value("optional", 0), Some(CifValue::Unknown)));
    assert!(matches!(test.value("late", 0), Some(CifValue::Unknown)));
    assert_eq!(test.text("late", 1), Some("yes"));
}

#[test]
fn writer_round_trips_categories_rows_and_xml_escaping() {
    let expected = parse_pdbml_document(XML.as_bytes()).expect("PDBML should parse");
    let xml = write_pdbml(&expected).expect("document should write");
    let actual = parse_pdbml_document(xml.as_bytes()).expect("written PDBML should parse");
    let expected_block = expected.first_block().expect("block should exist");
    let actual_block = actual.first_block().expect("block should exist");
    assert_eq!(actual_block.name(), expected_block.name());
    assert_eq!(
        actual_block
            .category("test")
            .map(crate::Category::row_count),
        Some(2)
    );
    assert_eq!(
        actual_block
            .category("entry")
            .and_then(|value| value.text("details", 0)),
        Some("A & B")
    );
}

#[test]
fn official_wwpdb_fixture_lowers_when_configured() {
    let Some(path) = std::env::var_os("PDBIOX_PDBML_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).expect("configured wwPDB fixture should be readable");
    let (_, structure, _) =
        read_pdbml(&bytes, &ReadOptions::new()).expect("wwPDB PDBML should lower");
    assert!(structure.atom_count() > 100);
    assert_eq!(structure.data().entry.id.as_deref(), Some("1CRN"));
}
