use super::*;
use crate::read;
use molframe_core::io::{InputBuffer, ReadOptions};

fn parse(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

#[test]
fn a_file_the_format_holds_exactly_reports_no_loss() {
    let source = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
END
";
    assert!(projection_loss(&parse(source), &PdbOptions::new()).is_empty());
}

#[test]
fn a_value_is_rounded_exactly_when_the_columns_cannot_hold_it() {
    assert!(rounded(1.2345, 1e3));
    assert!(!rounded(1.234, 1e3));
    assert!(rounded(0.125, 1e2));
    assert!(!rounded(0.12, 1e2));
}

const ATOM: &str = "ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N";

fn mentions(findings: &[Diagnostic], needle: &str) -> bool {
    findings.iter().any(|f| f.message().contains(needle))
}

#[test]
fn a_header_the_columns_hold_exactly_reports_no_loss() {
    let source = format!(
        "REMARK   2 RESOLUTION.    1.80 ANGSTROMS.\n\
         CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1\n{ATOM}1+\nTER       2      GLY A   1\nEND\n"
    );
    assert!(projection_loss(&parse(&source), &PdbOptions::new()).is_empty());
}

#[test]
fn a_resolution_and_cell_are_reported_only_when_rounded() {
    let source = format!(
        "REMARK   2 RESOLUTION.    1.80 ANGSTROMS.\n\
         CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1\n{ATOM}\nEND\n"
    );
    let mut structure = parse(&source);
    let findings = projection_loss(&structure, &PdbOptions::new());
    assert!(!mentions(&findings, "resolution"));
    assert!(!mentions(&findings, "cell"));
    let mut data = structure.data().clone();
    data.entry.resolution = Some(1.8049);
    if let Some(cell) = data.cell.as_mut() {
        cell.lengths[0] = 10.0004;
    }
    structure = Structure::new(data);
    let findings = projection_loss(&structure, &PdbOptions::new());
    assert!(mentions(&findings, "resolution"));
    assert!(mentions(&findings, "cell"));
}

#[test]
fn a_serial_is_reported_only_when_the_writer_renumbers_it() {
    let kept = projection_loss(&parse(&format!("{ATOM}\nEND\n")), &PdbOptions::new());
    assert!(!mentions(&kept, "serial"));
    let moved = ATOM.replacen("ATOM      1", "ATOM      7", 1);
    let findings = projection_loss(&parse(&format!("{moved}\nEND\n")), &PdbOptions::new());
    assert!(mentions(&findings, "7 becomes 1"));
}

#[test]
fn an_explicit_zero_charge_is_reported_and_a_signed_one_is_not() {
    let zero = projection_loss(&parse(&format!("{ATOM}0+\nEND\n")), &PdbOptions::new());
    assert!(mentions(&zero, "explicit zero charge"));
    let signed = projection_loss(&parse(&format!("{ATOM}1+\nEND\n")), &PdbOptions::new());
    assert!(!mentions(&signed, "zero charge"));
    let blank = projection_loss(&parse(&format!("{ATOM}\nEND\n")), &PdbOptions::new());
    assert!(!mentions(&blank, "zero charge"));
}
