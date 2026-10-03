/// Exact decimal fields parse to these values within rounding.
fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (left, right) in actual.iter().zip(expected) {
        assert!((left - right).abs() < 1e-9, "{actual:?} vs {expected:?}");
    }
}

use crate::{PdbHeaders, PdbHeadersExt};
use molframe_bench::{Sample, input};
use molframe_core::io::ReadOptions;

fn headers_of_text(text: &str) -> PdbHeaders {
    let body =
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 10.00           C\nEND\n";
    let source = format!("{text}{body}");
    let (structure, _) = crate::read(
        &molframe_core::io::InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("fixture reads");
    structure.pdb_headers().cloned().expect("headers kept")
}

fn headers_of(sample: Sample) -> PdbHeaders {
    let bytes = sample.pdb().expect("the sample ships a PDB");
    let (structure, _) = crate::read(&input(bytes), &ReadOptions::new()).expect("sample reads");
    structure.pdb_headers().cloned().expect("headers kept")
}

#[test]
fn resolution_is_read_when_deposited_and_absent_for_nmr() {
    assert_eq!(headers_of(Sample::Medium).resolution(), Some(1.74));
    assert_eq!(headers_of(Sample::Ensemble).resolution(), None);
}

#[test]
fn hemoglobin_declares_one_tetrameric_biomolecule() {
    let assemblies = headers_of(Sample::Medium).biomolecules();
    assert_eq!(assemblies.len(), 1);
    assert_eq!(assemblies[0].id, 1);
    assert_eq!(assemblies[0].groups.len(), 1);
    let group = &assemblies[0].groups[0];
    let chains: Vec<_> = group.chains.iter().map(|chain| &**chain).collect();
    assert_eq!(chains, ["A", "B", "C", "D"]);
    assert_eq!(group.operations.len(), 1);
    assert_close(&group.operations[0].rotation[0], &[1.0, 0.0, 0.0]);
    assert_close(&group.operations[0].translation, &[0.0; 3]);
}

#[test]
fn long_chain_lists_and_several_operators_group_correctly() {
    let text = "\
REMARK 350 BIOMOLECULE: 1                                             \n\
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C                     \n\
REMARK 350                    AND CHAINS: D, E                        \n\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000  \n\
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       10.00000  \n\
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000       20.00000  \n\
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000       30.00000  \n\
REMARK 350 BIOMOLECULE: 2                                             \n\
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A                           \n\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000  \n";
    let assemblies = headers_of_text(text).biomolecules();
    assert_eq!(assemblies.len(), 2);
    let first = &assemblies[0].groups[0];
    assert_eq!(first.chains.len(), 5);
    assert_eq!(first.operations.len(), 2);
    assert_eq!(first.operations[1].serial, 2);
    assert_close(&first.operations[1].rotation[1], &[0.0, -1.0, 0.0]);
    assert_close(&first.operations[1].translation, &[10.0, 20.0, 30.0]);
    assert_eq!(assemblies[1].groups[0].operations.len(), 1);
}

#[test]
fn missing_residues_are_read_from_the_remark_table_only() {
    let text = "\
REMARK 465                                                            \n\
REMARK 465 MISSING RESIDUES                                           \n\
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE             \n\
REMARK 465   M RES C SSSEQI                                           \n\
REMARK 465     MET A     1                                            \n\
REMARK 465     GLY A    12A                                           \n";
    let missing = headers_of_text(text).missing_residues();
    assert_eq!(missing.len(), 2);
    assert_eq!(&*missing[0].residue.name, "MET");
    assert_eq!(missing[0].residue.sequence, 1);
    assert_eq!(missing[0].model, None);
    assert_eq!(missing[1].residue.sequence, 12);
    assert_eq!(missing[1].residue.insertion, Some('A'));
}

#[test]
fn the_reader_fills_the_entry_resolution_from_remark_2() {
    let bytes = Sample::Medium.pdb().expect("4HHB ships a PDB");
    let (structure, _) = crate::read(&input(bytes), &ReadOptions::new()).expect("reads");
    assert_eq!(structure.data().entry.resolution, Some(1.74));
    let bytes = Sample::Ensemble.pdb().expect("2M7C ships a PDB");
    let (structure, _) = crate::read(&input(bytes), &ReadOptions::new()).expect("reads");
    assert_eq!(structure.data().entry.resolution, None);
}
