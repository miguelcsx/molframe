use crate::read;
use crate::writer::{PdbOptions, write};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const ATOMS: &str = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
END
";

fn parse(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

fn written(structure: &Structure) -> String {
    match write(structure, &PdbOptions::new()) {
        Ok(out) => out,
        Err(refusals) => panic!("write refused: {refusals:?}"),
    }
}

#[test]
fn a_title_read_from_continuation_lines_is_rewrapped_and_reads_back_whole() {
    let source = format!(
        "TITLE     {}\nTITLE    2 {}\nTITLE    3 {}\n{ATOMS}",
        "alpha ".repeat(11).trim_end(),
        "bravo ".repeat(11).trim_end(),
        "charlie delta echo",
    );
    let first = parse(&source);
    let title = first.data().entry.title.clone();
    let out = written(&first);
    assert!(out.lines().all(|line| line.len() <= 80));
    assert!(out.contains("TITLE    2 "));
    assert_eq!(parse(&out).data().entry.title, title);
}

#[test]
fn formal_charges_in_columns_79_and_80_survive_a_round_trip() {
    let source = "\
HETATM    1 FE    FE A   1      27.340  24.430   2.614  1.00 10.00          FE2+
HETATM    2  O   HOH A   2      26.266  25.413   2.842  1.00 11.00           O1-
HETATM    3  C   LIG A   3      26.335  27.770   3.258  1.00 14.00           C
END
";
    let out = written(&parse(source));
    let again = parse(&out);
    let charges: Vec<Option<i8>> = (0..3)
        .map(|n| {
            again
                .atom(molframe_core::index::AtomIndex::new(n))
                .and_then(molframe_core::structure::AtomRef::formal_charge)
        })
        .collect();
    assert_eq!(charges, [Some(2), Some(-1), None]);
    assert!(out.contains("2+\n"));
}
