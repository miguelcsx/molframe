use super::*;
use super::{DIPEPTIDE, parse};
use molframe_core::index::ModelIndex;
use molframe_core::structure::{AtomRef, ChainRef, Structure};

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

fn two_chain_model(number: u32, y: &str) -> String {
    use std::fmt::Write as _;
    let rows = [
        ("A", 1, "0.000"),
        ("A", 2, "3.800"),
        ("B", 1, "20.000"),
        ("B", 2, "23.800"),
    ];
    let mut text = format!("MODEL     {number:>4}\n");
    for (index, (chain, residue, x)) in rows.iter().enumerate() {
        let serial = index + 1;
        let _ = writeln!(
            text,
            "ATOM  {serial:>5}  CA  ALA {chain}{residue:>4}    {x:>8}{y:>8}   0.000  1.00  0.00           C"
        );
        if *residue == 2 {
            text.push_str("TER\n");
        }
    }
    text.push_str("ENDMDL\n");
    text
}

/// Later models only add coordinates. Their `TER` and `ENDMDL` records must not
/// touch the residue ranges the first model fixed, which used to fail with a
/// reversed-range consistency error as soon as a file had two chains.
#[test]
fn multi_model_files_with_several_chains_keep_the_first_model_hierarchy() {
    let text = format!(
        "{}{}{}END\n",
        two_chain_model(1, "1.000"),
        two_chain_model(2, "2.000"),
        two_chain_model(3, "3.000")
    );
    let (structure, _) = parse(&text);
    assert_eq!(structure.model_count(), 3);
    assert_eq!(structure.chain_count(), 2);
    assert_eq!(structure.residue_count(), 4);
    assert_eq!(structure.atom_count(), 4);
    let residues: Vec<_> = (0..4)
        .map(|index| {
            structure
                .data()
                .topology
                .residues
                .atoms(molframe_core::index::ResidueIndex::new(index))
        })
        .collect();
    assert_eq!(
        residues,
        vec![Some(0..1), Some(1..2), Some(2..3), Some(3..4)]
    );
}
