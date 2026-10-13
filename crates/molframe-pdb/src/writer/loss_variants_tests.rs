use super::*;
use molframe_core::{Format, InputBuffer, ReadOptions};

const PQR: &str = "ATOM      1  N   ALA A   1      11.104  13.207   9.124 -0.3000 1.5500\nEND\n";
const PDBQT: &str = "ROOT\n\
ATOM      1  N   LIG A   1      11.104  13.207   9.124  1.00  0.00    -0.300 N \n\
ENDROOT\nTORSDOF 0\n";

fn read(text: &str, format: Format) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().format(format);
    let result = if format == Format::Pqr {
        crate::read_pqr(&input, &options)
    } else {
        crate::read_pdbqt(&input, &options)
    };
    match result {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

fn messages(findings: &[Diagnostic]) -> String {
    findings
        .iter()
        .map(|d| d.message().to_owned())
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn canonical_input_is_reported_lossless() {
    let pqr = pqr_loss(&read(PQR, Format::Pqr), &PdbOptions::new());
    assert!(pqr.is_empty(), "{}", messages(&pqr));
    let pdbqt = pdbqt_loss(&read(PDBQT, Format::Pdbqt), &PdbOptions::new());
    assert!(pdbqt.is_empty(), "{}", messages(&pdbqt));
}

#[test]
fn pqr_reports_a_charge_and_radius_beyond_four_decimals() {
    let text = PQR.replace("-0.3000 1.5500", "-0.30004 1.55004");
    let findings = pqr_loss(&read(&text, Format::Pqr), &PdbOptions::new());
    let text = messages(&findings);
    assert!(text.contains("partial charge"), "{text}");
    assert!(text.contains("radius"), "{text}");
}

#[test]
fn pdbqt_reports_a_charge_beyond_three_decimals() {
    let text = PDBQT.replace("-0.300 N", "-0.3004 N");
    let findings = pdbqt_loss(&read(&text, Format::Pdbqt), &PdbOptions::new());
    assert!(messages(&findings).contains("partial charge"));
}

#[test]
fn pdbqt_reports_a_type_wider_than_its_columns() {
    let text = PDBQT.replace("-0.300 N ", "-0.300 NAA");
    let findings = pdbqt_loss(&read(&text, Format::Pdbqt), &PdbOptions::new());
    assert!(messages(&findings).contains("AutoDock type"));
}

#[test]
fn pqr_reports_the_occupancy_it_has_no_column_for() {
    let structure = read(PDBQT, Format::Pdbqt);
    let findings = pqr_loss(&structure, &PdbOptions::new());
    let text = messages(&findings);
    assert!(!text.contains("occupancy") || text.contains("PQR has no column"));
}

#[test]
fn a_torsion_tree_is_reported_on_read() {
    let tree = "ROOT\n\
ATOM      1  C   LIG A   1      11.104  13.207   9.124  1.00  0.00    +0.100 C \n\
ENDROOT\nBRANCH   1   2\n\
ATOM      2  C   LIG A   1      12.104  13.207   9.124  1.00  0.00    +0.100 C \n\
ENDBRANCH   1   2\nTORSDOF 1\n";
    let input = InputBuffer::from_bytes(tree.as_bytes().to_vec());
    let options = ReadOptions::new().format(Format::Pdbqt);
    let Ok((_, findings)) = crate::read_pdbqt(&input, &options) else {
        panic!("read failed")
    };
    assert!(messages(&findings).contains("torsion tree"));
    let (_, quiet) = match crate::read_pdbqt(
        &InputBuffer::from_bytes(PDBQT.as_bytes().to_vec()),
        &options,
    ) {
        Ok(read) => read,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    assert!(!messages(&quiet).contains("torsion tree"));
}

#[test]
fn the_pqr_reader_does_not_mistake_charge_columns_for_occupancy() {
    let structure = read(PQR, Format::Pqr);
    let Some(atom) = structure.data().atoms().next() else {
        panic!("no atom")
    };
    assert_eq!(atom.occupancy(), None);
    assert_eq!(atom.b_factor(), None);
}

#[test]
fn pdbqt_charge_and_type_land_in_the_published_columns() {
    let written = match crate::write_pdbqt(&read(PDBQT, Format::Pdbqt), &PdbOptions::new()) {
        Ok(written) => written,
        Err(findings) => panic!("write failed: {findings:?}"),
    };
    let Some(row) = written.lines().find(|line| line.starts_with("ATOM")) else {
        panic!("no atom row")
    };
    assert_eq!(row.get(70..76), Some("-0.300"));
    assert_eq!(row.get(77..79), Some("N "));
}
