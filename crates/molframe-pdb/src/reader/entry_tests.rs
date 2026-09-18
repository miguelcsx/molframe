use super::*;
use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef};

fn parse(text: &str) -> (Structure, Vec<Diagnostic>) {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

const DIPEPTIDE: &str = "\
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

#[test]
fn several_models_share_one_topology_and_contribute_frames() {
    let text = "\
MODEL        1
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  N   GLY A   1      1.500   1.500   1.500  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.500   2.500   2.500  1.00  0.00           C
ENDMDL
END
";
    let (structure, _) = parse(text);
    assert_eq!(structure.model_count(), 2);
    assert_eq!(structure.atom_count(), 2, "the topology is not duplicated");

    let first = structure
        .model_positions(ModelIndex::new(0))
        .map(<[[f32; 3]]>::len);
    let second = structure
        .model_positions(ModelIndex::new(1))
        .map(<[[f32; 3]]>::len);
    assert_eq!((first, second), (Some(2), Some(2)));
    let numbers: Vec<_> = structure
        .data()
        .models()
        .filter_map(molframe_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [1, 2]);
}

#[test]
fn reading_only_the_first_pdb_model_does_not_append_an_empty_frame() {
    let text = "\
MODEL        4
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ENDMDL
MODEL        8
ATOM      1  N   GLY A   1      2.000   2.000   2.000  1.00  0.00           N
ENDMDL
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().only_first_model(true);
    let (structure, _) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    assert_eq!(structure.model_count(), 1);
    let number = structure
        .data()
        .models()
        .next()
        .and_then(molframe_core::structure::ModelRef::number);
    assert_eq!(number, Some(4));
}

#[test]
fn a_later_pdb_model_with_different_atom_identity_is_ragged() {
    let text = "\
MODEL        1
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  N   GLY A   1      1.500   1.500   1.500  1.00  0.00           N
ATOM      2  O   GLY A   1      2.500   2.500   2.500  1.00  0.00           O
ENDMDL
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().mode(molframe_core::io::ParseMode::Recover);
    let (structure, _) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("ragged read failed: {findings:?}"),
    };
    let Some(models) = structure.ragged_models() else {
        panic!("identity-changing models must be ragged")
    };
    let names: Vec<Vec<_>> = models
        .iter()
        .map(|model| model.data().atoms().filter_map(AtomRef::name).collect())
        .collect();
    assert_eq!(names, [["N", "CA"], ["N", "O"]]);
    let numbers: Vec<_> = structure
        .data()
        .models()
        .filter_map(molframe_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [1, 2]);
}

#[test]
fn pdb_models_with_different_atom_counts_are_ragged() {
    let text = "MODEL        1\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N\nATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C\nENDMDL\nMODEL        2\nATOM      1  N   GLY A   1      1.500   1.500   1.500  1.00  0.00           N\nENDMDL\nEND\n";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (structure, _) = match read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("ragged read failed: {findings:?}"),
    };
    let counts: Vec<_> = structure
        .ragged_models()
        .into_iter()
        .flatten()
        .map(Structure::atom_count)
        .collect();
    assert_eq!(counts, [2, 1]);
}

#[test]
fn a_pdb_altloc_component_identity_is_kept_per_atom() {
    let text = "\
ATOM      1  N  AGLY A   1      1.000   1.000   1.000  0.50  0.00           N
ATOM      2  CA BALA A   1      2.000   2.000   2.000  0.50  0.00           C
END
";
    let (structure, findings) = parse(text);
    let components: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(AtomRef::component_name)
        .collect();
    assert_eq!(components, ["GLY", "ALA"]);
    assert!(findings.iter().any(|finding| finding.code() == Code::W3012));
}

#[test]
fn a_missing_element_column_is_inferred_and_the_inference_is_reported() {
    let text = "\
ATOM      1  CA  GLY A   1      1.000   1.000   1.000  1.00  0.00
HETATM    2 ZN    ZN A 100      2.000   2.000   2.000  1.00  0.00
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new()
        .missing_element_policy(molframe_core::MissingElementPolicy::InferFromAtomName);
    let (structure, findings) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let elements: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(AtomRef::element)
        .map(molframe_core::Element::symbol)
        .collect();
    assert_eq!(
        elements,
        ["C", "Zn"],
        "the leading space distinguishes the two"
    );
    assert!(findings.iter().all(|finding| finding.code() == Code::W3203));
    assert_eq!(findings.len(), 2);
}

#[test]
fn a_missing_element_is_unknown_without_an_explicit_inference_policy() {
    let text = "\
ATOM      1  CA  GLY A   1      1.000   1.000   1.000  1.00  0.00
END
";
    let (structure, findings) = parse(text);
    assert_eq!(
        structure.data().atoms().next().and_then(AtomRef::element),
        Some(molframe_core::Element::UNKNOWN)
    );
    assert!(findings.iter().any(|finding| finding.code() == Code::W3203));
}

#[test]
fn a_serial_past_the_decimal_range_reads_through_the_extended_scheme() {
    let text = "\
ATOM  A0000  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
END
";
    let (structure, _) = parse(text);
    let chunk = structure
        .data()
        .chunks
        .first()
        .map(|chunk| chunk.atom_site_id(0));
    assert_eq!(chunk, Some(Some(100_000)));
}

#[test]
fn an_unreadable_coordinate_is_reported_and_the_rest_of_the_file_still_reads() {
    let text = "\
ATOM      1  N   GLY A   1       abcd   1.000   1.000  1.00  0.00           N
ATOM      2  CA  GLY A   1      2.000   2.000   2.000  1.00  0.00           C
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().mode(molframe_core::io::ParseMode::Recover);
    let (structure, findings) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("recover mode refused: {findings:?}"),
    };
    assert_eq!(
        structure.atom_count(),
        2,
        "the row survives without its position"
    );
    assert!(findings.iter().any(|finding| finding.code() == Code::E1202));

    let positions: Vec<_> = structure.data().atoms().map(AtomRef::position).collect();
    assert_eq!(positions.first().copied(), Some(None));
    assert!(positions.get(1).copied().flatten().is_some());
}

#[test]
fn an_invalidating_coordinate_refuses_permissive_mode() {
    let text = "\
ATOM      1  N   GLY A   1       abcd   1.000   1.000  1.00  0.00           N
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let refused = read(&input, &ReadOptions::new());
    assert!(refused.is_err());
}

#[test]
fn a_file_with_no_coordinate_records_is_refused_rather_than_read_as_empty() {
    let input = InputBuffer::from_bytes(b"HEADER    NOTHING\nEND\n".to_vec());
    let refused = read(&input, &ReadOptions::new());
    let codes: Vec<_> = refused
        .err()
        .unwrap_or_default()
        .iter()
        .map(Diagnostic::code)
        .collect();
    assert!(codes.contains(&Code::E1001));
}

#[test]
fn hydrogens_are_dropped_when_the_caller_asks_for_that() {
    let text = "\
ATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N
ATOM      2  H   GLY A   1      2.000   2.000   2.000  1.00  0.00           H
END
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().discard_hydrogens(true);
    let read = read(&input, &options).map(|(structure, _)| structure.atom_count());
    assert_eq!(read.ok(), Some(1));
}

#[test]
fn every_residue_covers_the_atoms_that_were_read_into_it() {
    let (structure, _) = parse(DIPEPTIDE);
    assert!(
        molframe_core::structure::validate(structure.data()).is_empty(),
        "the structure a read produces must satisfy its own invariants"
    );

    let counts: Vec<_> = structure
        .data()
        .residues()
        .map(|residue| residue.atoms().count())
        .collect();
    assert_eq!(counts, [4, 2]);
}

#[test]
fn walking_the_hierarchy_reaches_a_named_atom_of_a_named_residue() {
    let (structure, _) = parse(DIPEPTIDE);
    let atom = structure
        .data()
        .model(ModelIndex::new(0))
        .and_then(|model| model.chain("A"))
        .and_then(|chain| chain.residue(2))
        .and_then(|residue| residue.atom("CA"));
    assert!(atom.is_some());
    assert_eq!(
        structure
            .data()
            .chains()
            .next()
            .map(|chain: ChainRef<'_>| chain.residues().count()),
        Some(2)
    );
}
