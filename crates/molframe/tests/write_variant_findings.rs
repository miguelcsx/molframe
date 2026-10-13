//! The PQR and PDBQT writers report what they round or drop, and stay silent
//! on rows they hold exactly.

use molframe::{ReadOptions, Structure};
use molframe_pdb::PdbOptions;

const PQR: &str = "ATOM      1  N   ALA A   1      11.104  13.207   9.124 -0.3000 1.5500\nEND\n";
const PDBQT: &str = "ROOT\n\
ATOM      1  N   LIG A   1      11.104  13.207   9.124  1.00  0.00    -0.300 N \n\
ENDROOT\nTORSDOF 0\n";

fn read(text: &str, name: &str) -> Structure {
    match molframe::read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(failure) => panic!("{name}: {failure:?}"),
    }
}

fn codes(findings: &molframe::Findings) -> Vec<String> {
    findings
        .to_vec()
        .iter()
        .map(|finding| format!("{} {}", finding.code(), finding.message()))
        .collect()
}

#[test]
fn a_canonical_pqr_writes_without_findings() {
    let structure = read(PQR, "in.pqr");
    let Ok((text, findings)) = molframe::write_pqr_with_findings(&structure, &PdbOptions::new())
    else {
        panic!("pqr write failed");
    };
    assert_eq!(text, PQR.to_owned());
    assert!(codes(&findings).is_empty(), "{:?}", codes(&findings));
}

#[test]
fn a_pqr_charge_beyond_its_decimals_is_reported() {
    let structure = read(&PQR.replace("-0.3000", "-0.30004"), "in.pqr");
    let Ok((_, findings)) = molframe::write_pqr_with_findings(&structure, &PdbOptions::new())
    else {
        panic!("pqr write failed");
    };
    assert!(
        codes(&findings)
            .iter()
            .any(|m| m.contains("W4101") && m.contains("charge"))
    );
}

#[test]
fn a_canonical_pdbqt_writes_without_findings() {
    let structure = read(PDBQT, "in.pdbqt");
    let Ok((text, findings)) = molframe::write_pdbqt_with_findings(&structure, &PdbOptions::new())
    else {
        panic!("pdbqt write failed");
    };
    assert_eq!(text, PDBQT.to_owned());
    assert!(codes(&findings).is_empty(), "{:?}", codes(&findings));
}
