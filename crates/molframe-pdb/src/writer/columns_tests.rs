use super::atom_name_field;
use crate::read;
use crate::writer::{PdbOptions, write};
use molframe_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
TER       2      GLY A   1
HETATM    3 FE   HEM A 142      18.362  18.488  23.755  1.00 18.07          FE
HETATM    4  O   HOH A 201      26.266  25.413   2.842  1.00 11.00           O
END
";

fn rewritten() -> String {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let Ok((structure, _)) = read(&input, &ReadOptions::new()) else {
        panic!("read failed");
    };
    match write(&structure, &PdbOptions::new()) {
        Ok(out) => out,
        Err(refusals) => panic!("write refused: {refusals:?}"),
    }
}

#[test]
fn a_two_letter_element_starts_its_name_in_the_first_name_column() {
    assert_eq!(atom_name_field("FE", "FE"), "FE  ");
    assert_eq!(atom_name_field("CA", "C"), " CA ");
    assert_eq!(atom_name_field("CA", "CA"), "CA  ");
}

#[test]
fn metals_keep_their_upper_case_element_and_name_columns() {
    let out = rewritten();
    assert!(out.contains("HETATM    3 FE   HEM A 142"));
    assert!(
        out.lines()
            .any(|line| line.starts_with("HETATM    3") && line.ends_with("FE"))
    );
}

#[test]
fn only_a_polymer_chain_is_terminated_and_the_record_names_its_last_residue() {
    let out = rewritten();
    let terminators: Vec<&str> = out.lines().filter(|l| l.starts_with("TER")).collect();
    assert_eq!(terminators, ["TER       2      GLY A   1"]);
}
