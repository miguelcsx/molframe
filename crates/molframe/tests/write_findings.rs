//! The CIF-family and PDB writers report what they rewrite or cannot carry,
//! and stay silent on input the target holds exactly.

use molframe::{ReadOptions, Structure};
use molframe_cif::CifWriteOptions;
use molframe_pdb::PdbOptions;

const POLYMER: &str = "\
ATOM      1  N   ALA A  10      11.000  12.000  13.000  1.00 10.00           N
ATOM      2  CA  ALA A  10      12.000  12.000  13.000  1.00 10.00           C
TER       3      ALA A  10
";
const WATER: &str =
    "HETATM    4  O   HOH A 101      20.000  20.000  20.000  1.00 20.00           O\n";
const INSERTED: &str = "\
ATOM      3  N   GLY A  10A     14.000  12.000  13.000  1.00 10.00           N
";

fn read(text: &str, name: &str) -> Structure {
    match molframe::read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(failure) => panic!("{name}: {failure:?}"),
    }
}

fn messages(findings: &molframe::Findings) -> Vec<String> {
    findings
        .to_vec()
        .iter()
        .map(|finding| format!("{} {}", finding.code(), finding.message()))
        .collect()
}

fn cif_of(pdb: &str) -> (String, Vec<u8>, Vec<String>, Vec<String>) {
    let structure = read(pdb, "in.pdb");
    let options = CifWriteOptions::new();
    let Ok((text, text_findings)) = molframe::write_mmcif_with_findings(&structure, &options)
    else {
        panic!("mmcif write failed");
    };
    let Ok((bytes, binary_findings)) = molframe::write_bcif_with_findings(&structure, &options)
    else {
        panic!("bcif write failed");
    };
    (
        text,
        bytes,
        messages(&text_findings),
        messages(&binary_findings),
    )
}

fn any(findings: &[String], needle: &str) -> bool {
    findings.iter().any(|message| message.contains(needle))
}

#[test]
fn a_repeated_chain_label_is_reported_by_both_cif_writers() {
    let (_, _, text, binary) = cif_of(&format!("{POLYMER}{WATER}END\n"));
    for findings in [&text, &binary] {
        assert!(any(findings, "W4102"), "{findings:?}");
        assert!(any(findings, "`A` becomes `A2`"), "{findings:?}");
    }
}

#[test]
fn a_shifted_residue_number_is_reported_by_both_cif_writers() {
    let polymer = POLYMER.replacen("TER", &format!("{INSERTED}TER"), 1);
    let (_, _, text, binary) = cif_of(&format!("{polymer}END\n"));
    for findings in [&text, &binary] {
        assert!(
            any(findings, "shifted normalised sequence number"),
            "{findings:?}"
        );
    }
}

#[test]
fn distinct_identifiers_are_not_reported_as_rewritten() {
    let (_, _, text, binary) = cif_of(&format!("{POLYMER}END\n"));
    for findings in [&text, &binary] {
        assert!(!any(findings, "W4102"), "{findings:?}");
    }
}

#[test]
fn a_cif_read_back_from_cif_writes_with_no_findings() {
    let (text, _, _, _) = cif_of(&format!("{POLYMER}END\n"));
    let again = read(&text, "again.cif");
    let options = CifWriteOptions::new();
    let Ok((_, text_findings)) = molframe::write_mmcif_with_findings(&again, &options) else {
        panic!("mmcif write failed");
    };
    let Ok((_, binary_findings)) = molframe::write_bcif_with_findings(&again, &options) else {
        panic!("bcif write failed");
    };
    assert!(
        text_findings.to_vec().is_empty(),
        "{:?}",
        messages(&text_findings)
    );
    assert!(
        binary_findings.to_vec().is_empty(),
        "{:?}",
        messages(&binary_findings)
    );
}

#[test]
fn the_plain_writers_agree_with_the_ones_that_report() {
    let structure = read(&format!("{POLYMER}{WATER}END\n"), "in.pdb");
    let Ok(plain) = molframe::write_mmcif(&structure) else {
        panic!("plain write failed");
    };
    let Ok((reported, _)) =
        molframe::write_mmcif_with_findings(&structure, &CifWriteOptions::new())
    else {
        panic!("reporting write failed");
    };
    assert_eq!(plain, reported);
}

#[test]
fn pdb_only_extensions_are_reported_when_written_to_cif() {
    let pdb = format!("TITLE     A TEST\nREMARK   3 REFINEMENT.\n{POLYMER}END\n");
    let structure = read(&pdb, "in.pdb");
    let has_headers = structure
        .engine()
        .extensions()
        .keys()
        .any(|key| key == molframe_pdb::PDB_HEADERS_EXTENSION);
    let Ok((_, findings)) =
        molframe::write_mmcif_with_findings(&structure, &CifWriteOptions::new())
    else {
        panic!("mmcif write failed");
    };
    assert!(has_headers);
    assert!(any(
        &messages(&findings),
        molframe_pdb::PDB_HEADERS_EXTENSION
    ));
}

#[test]
fn a_pdb_write_reports_each_header_and_serial_loss_once() {
    let source = format!(
        "REMARK   2 RESOLUTION.    1.80 ANGSTROMS.\n{}END\n",
        POLYMER.replacen("ATOM      1", "ATOM      5", 1)
    );
    let structure = read(&source, "in.pdb");
    let Ok((_, findings)) = molframe::write_pdb_with_findings(&structure, &PdbOptions::new())
    else {
        panic!("pdb write failed");
    };
    let findings = messages(&findings);
    assert!(any(&findings, "5 becomes 1"), "{findings:?}");
    assert!(!any(&findings, "resolution"), "{findings:?}");
}
