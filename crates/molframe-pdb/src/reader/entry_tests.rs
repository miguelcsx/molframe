use super::*;
use molframe_core::structure::Structure;
use molframe_core::structure::{AtomRef, ResidueRef};

pub(super) fn parse(text: &str) -> (Structure, Vec<Diagnostic>) {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

pub(super) const DIPEPTIDE: &str = "\
HEADER    HYDROLASE                               01-JAN-00   1ABC
CRYST1   61.500   61.500   97.300  90.00  90.00 120.00 P 32 2 1
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      26.266  25.413   2.842  1.00 11.00           C
ATOM      3  C   GLY A   1      26.913  26.639   3.531  1.00 12.00           C
ATOM      4  O   GLY A   1      27.886  26.463   4.263  1.00 13.00           O
ATOM      5  N   ASN A   2      26.335  27.770   3.258  1.00 14.00           N
ATOM      6  CA  ASN A   2      26.850  29.021   3.898  1.00 15.00           C
TER       7      ASN A   2
END
";

#[test]
fn a_dipeptide_reads_with_the_hierarchy_the_file_describes() {
    let (structure, _) = parse(DIPEPTIDE);
    assert_eq!(structure.atom_count(), 6);
    assert_eq!(structure.residue_count(), 2);
    assert_eq!(structure.chain_count(), 1);
    assert_eq!(structure.model_count(), 1);
}

#[test]
fn atom_records_do_not_invent_polymer_chemistry() {
    let (structure, _) = parse(DIPEPTIDE);
    let Some(chain) = structure.data().chains().next() else {
        panic!("chain absent");
    };

    assert_eq!(
        chain.polymer_kind(),
        molframe_core::topology::PolymerKind::None
    );
}

#[test]
fn the_entry_identifier_and_cell_are_read_from_their_own_records() {
    let (structure, _) = parse(DIPEPTIDE);
    assert_eq!(structure.data().entry.id.as_deref(), Some("1ABC"));
    let cell = structure.data().cell;
    assert!(cell.is_some_and(|cell| (cell.lengths[0] - 61.5).abs() < 1e-6));
    assert!(cell.is_some_and(|cell| (cell.angles[2] - 120.0).abs() < 1e-6));
}

#[test]
fn every_format_33_metadata_record_is_preserved_in_order() {
    use crate::{PdbHeadersExt, PdbOptions, write};
    use std::fmt::Write as _;

    let names = [
        "HEADER", "OBSLTE", "TITLE", "SPLIT", "CAVEAT", "COMPND", "SOURCE", "KEYWDS", "EXPDTA",
        "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK", "DBREF", "DBREF1",
        "DBREF2", "SEQADV", "SEQRES", "MODRES", "HET", "HETNAM", "HETSYN", "FORMUL", "HELIX",
        "SHEET", "SSBOND", "LINK", "CISPEP", "SITE", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
        "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3",
    ];
    let mut source = String::new();
    for name in names {
        match name {
            "HEADER" => source
                .push_str("HEADER    TEST CLASSIFICATION                     01-JAN-00   1ABC\n"),
            "CRYST1" => {
                source.push_str("CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1\n");
            }
            _ => {
                let _ = writeln!(source, "{name:<6}    deposited value");
            }
        }
    }
    source.push_str(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 10.00           C\nEND\n",
    );
    let (structure, _) = parse(&source);
    let Some(headers) = structure.pdb_headers() else {
        panic!("metadata extension absent");
    };
    let actual: Vec<_> = headers
        .records()
        .iter()
        .map(crate::PdbHeaderRecord::name)
        .collect();
    assert_eq!(actual, names);
    assert_eq!(headers.classification(), Some("TEST CLASSIFICATION"));
    assert_eq!(headers.deposition_date(), Some("01-JAN-00"));

    let written = write(&structure, &PdbOptions::new())
        .unwrap_or_else(|findings| panic!("write failed: {findings:?}"));
    let (round_trip, _) = parse(&written);
    let Some(round_trip_headers) = round_trip.pdb_headers() else {
        panic!("round-trip metadata extension absent");
    };
    let round_trip_names: Vec<_> = round_trip_headers
        .records()
        .iter()
        .map(crate::PdbHeaderRecord::name)
        .collect();
    assert_eq!(round_trip_names, names);
}

#[test]
fn atom_names_elements_and_positions_survive_the_read() {
    let (structure, _) = parse(DIPEPTIDE);
    let names: Vec<_> = structure.data().atoms().filter_map(AtomRef::name).collect();
    assert_eq!(names, ["N", "CA", "C", "O", "N", "CA"]);

    let elements: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(AtomRef::element)
        .map(molframe_core::Element::symbol)
        .collect();
    assert_eq!(elements, ["N", "C", "C", "O", "N", "C"]);

    let first = structure.data().atoms().next().and_then(AtomRef::position);
    assert!(first.is_some_and(|position| (position[0] - 27.340).abs() < 1e-3));
}

#[test]
fn residues_carry_the_numbers_the_depositor_gave_them() {
    let (structure, _) = parse(DIPEPTIDE);
    let numbers: Vec<_> = structure
        .data()
        .residues()
        .filter_map(ResidueRef::auth_seq_id)
        .collect();
    assert_eq!(numbers, [1, 2]);
    let names: Vec<_> = structure
        .data()
        .residues()
        .filter_map(ResidueRef::name)
        .collect();
    assert_eq!(names, ["GLY", "ASN"]);
}

#[test]
fn conect_records_become_deduplicated_file_provenance_edges() {
    let source = format!("{DIPEPTIDE}CONECT    1    2    3\nCONECT    2    1\n");
    let (structure, _) = parse(&source);
    let bonds: Vec<_> = structure.data().bonds.iter().collect();
    assert_eq!(bonds.len(), 2);
    assert_eq!((bonds[0].atom_a.get(), bonds[0].atom_b.get()), (0, 1));
    assert_eq!(bonds[0].provenance, molframe_core::BondProvenance::File);
    assert_eq!(bonds[0].order, molframe_core::BondOrder::Unknown);
}

#[test]
fn an_insertion_code_starts_a_new_residue_rather_than_extending_one() {
    let text = "\
ATOM      1  N   HIS L 163      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  N   HIS L 163A     2.000   2.000   2.000  1.00  0.00           N
ATOM      3  N   HIS L 163B     3.000   3.000   3.000  1.00  0.00           N
END
";
    let (structure, _) = parse(text);
    assert_eq!(
        structure.residue_count(),
        3,
        "163, 163A and 163B are three residues"
    );
    let codes: Vec<_> = structure
        .data()
        .residues()
        .map(|residue| residue.ins_code().unwrap_or("").to_owned())
        .collect();
    assert_eq!(codes, ["", "A", "B"]);
}

#[test]
fn a_change_of_chain_starts_a_new_chain() {
    let text = "\
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  N   GLY B   1      2.000   2.000   2.000  1.00  0.00           N
END
";
    let (structure, _) = parse(text);
    assert_eq!(structure.chain_count(), 2);
    let labels: Vec<_> = structure
        .data()
        .chains()
        .filter_map(|chain| chain.auth_asym_id().and_then(|s| structure.resolve(s)))
        .map(str::to_owned)
        .collect();
    assert_eq!(labels, ["A", "B"]);
}

#[path = "entry_models_tests.rs"]
mod models;

#[test]
fn a_file_spanning_several_decode_windows_keeps_row_order_and_findings() {
    let rows = WINDOW_LINES * 2 + 17;
    let text = (0..rows)
        .map(|index| {
            let x = if index == WINDOW_LINES + 3 {
                "bad"
            } else {
                "27.340"
            };
            let serial = index % 90_000 + 1;
            let residue = index % 9_000 + 1;
            format!(
                "ATOM  {serial:>5}  CA  GLY A{residue:>4}    {x:>8}  24.430   2.614  1.00 10.00           C\n"
            )
        })
        .collect::<Vec<_>>()
        .concat();
    let input = InputBuffer::from_bytes(text.into_bytes());
    let Err(findings) = read(&input, &ReadOptions::new()) else {
        panic!("a malformed coordinate must fail the read");
    };
    assert_eq!(findings.len(), 1);
    let line = findings
        .first()
        .and_then(molframe_core::Diagnostic::span)
        .map(|span| span.start.line);
    assert_eq!(line, u64::try_from(WINDOW_LINES + 4).ok());
}
