//! End-to-end small-molecule chemistry on the standard amino-acid dictionary.

mod common;
use common::{amino_acid_ccd, atom_records, haemoglobin_dimer, rows};
use std::collections::BTreeSet;

fn ccd() -> [String; 4] {
    [
        "--ccd".to_owned(),
        amino_acid_ccd().display().to_string(),
        "--ccd-version".to_owned(),
        "wwPDB-2026-10-03".to_owned(),
    ]
}

fn invoke(directory: &std::path::Path, head: &[&str], tail: &[&str]) -> Vec<Vec<String>> {
    let ccd = ccd();
    let mut arguments = head.to_vec();
    arguments.extend(ccd.iter().map(String::as_str));
    arguments.extend_from_slice(tail);
    rows(directory, &arguments)
}

#[test]
fn the_benzo_ring_of_tryptophan_matches_in_twelve_ring_orderings() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = invoke(
        directory.path(),
        &[
            "chem",
            "smarts",
            "TRP",
            "--component",
            "--pattern",
            "c1ccccc1",
        ],
        &[],
    );
    assert_eq!(table[0], ["match", "query_atom", "atom"]);
    let matches: BTreeSet<&str> = table[1..].iter().map(|row| row[0].as_str()).collect();
    assert_eq!(matches.len(), 12, "six rotations of two directions");
    let ring = ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"];
    assert!(table[1..].iter().all(|row| ring.contains(&row[2].as_str())));
}

#[test]
fn side_chain_amides_are_found_once_per_asparagine_and_glutamine() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let dimer = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let expected = std::fs::read_to_string(&dimer)
        .expect("dimer")
        .lines()
        .filter(|line| line.starts_with("ATOM") && line[12..16].trim() == "CA")
        .filter(|line| matches!(line[17..20].trim(), "ASN" | "GLN"))
        .count();
    let dimer = dimer.display().to_string();
    let table = invoke(
        directory.path(),
        &["chem", "smarts", &dimer, "--pattern", "C(=O)N"],
        &[],
    );
    let matches: BTreeSet<&str> = table[1..].iter().map(|row| row[0].as_str()).collect();
    assert!(expected > 0);
    assert_eq!(matches.len(), expected);
}

#[test]
fn a_malformed_pattern_is_refused_with_its_position() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let ccd = ccd();
    let mut arguments = vec!["chem", "smarts", "TRP", "--component", "--pattern", "C(("];
    arguments.extend(ccd.iter().map(String::as_str));
    let output = common::run(directory.path(), &arguments);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("byte"));
}

fn charges(directory: &std::path::Path, input: &str, passes: &str) -> Vec<f64> {
    let table = invoke(
        directory,
        &["chem", "peoe", input],
        &[
            "--passes",
            passes,
            "--initial-damping",
            "0.5",
            "--damping-factor",
            "0.5",
            "--minimum-difference",
            "0",
            "--profile",
            "gasteiger-marsili",
        ],
    );
    assert_eq!(table[0], ["atom", "charge"]);
    table[1..]
        .iter()
        .map(|row| row[1].parse().expect("charge"))
        .collect()
}

#[test]
fn peoe_assigns_every_atom_a_finite_charge_and_conserves_the_total() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let dimer = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let atoms = atom_records(&dimer);
    let input = dimer.display().to_string();
    let few = charges(directory.path(), &input, "3");
    let many = charges(directory.path(), &input, "12");
    assert_eq!(few.len(), atoms);
    assert!(many.iter().all(|charge| charge.is_finite()));
    // Charge moves between atoms; it is never created, whatever the passes.
    let (total_few, total_many): (f64, f64) = (few.iter().sum(), many.iter().sum());
    assert!(
        (total_few - total_many).abs() < 1e-6,
        "{total_few} vs {total_many}"
    );
    assert!(
        few.iter().zip(&many).any(|(a, b)| (a - b).abs() > 1e-4),
        "more passes move more charge"
    );
    assert_eq!(
        many,
        charges(directory.path(), &input, "12"),
        "deterministic"
    );
}

#[test]
fn peoe_refuses_hydrogen_names_the_dictionary_does_not_define() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = common::hydrogenated_peptide(directory.path())
        .display()
        .to_string();
    let ccd = ccd();
    let mut arguments = vec![
        "chem",
        "peoe",
        peptide.as_str(),
        "--passes",
        "12",
        "--initial-damping",
        "0.5",
        "--damping-factor",
        "0.5",
        "--minimum-difference",
        "0",
        "--profile",
        "gasteiger-marsili",
    ];
    arguments.extend(ccd.iter().map(String::as_str));
    let output = common::run(directory.path(), &arguments);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("UnknownAtom"));
}
