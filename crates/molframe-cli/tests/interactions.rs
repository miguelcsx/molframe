//! End-to-end interaction commands on real hydrogenated peptide chemistry.

mod common;
use common::{amino_acid_ccd, column, hydrogenated_peptide, rows};

fn ccd() -> Vec<String> {
    vec![
        "--ccd".into(),
        amino_acid_ccd().display().to_string(),
        "--ccd-version".into(),
        "wwPDB-2026-10-03".into(),
    ]
}

fn invoke(directory: &std::path::Path, command: &[&str], extra: &[&str]) -> Vec<Vec<String>> {
    let peptide = hydrogenated_peptide(directory);
    let ccd = ccd();
    let mut arguments = command.to_vec();
    let peptide = peptide.display().to_string();
    arguments.push(&peptide);
    arguments.extend(ccd.iter().map(String::as_str));
    arguments.extend_from_slice(extra);
    rows(directory, &arguments)
}

#[test]
fn hydrogen_bonds_respect_the_stated_geometry() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = invoke(
        directory.path(),
        &["interactions", "hbonds"],
        &[
            "--max-donor-acceptor-distance",
            "3.5",
            "--min-angle-degrees",
            "120",
        ],
    );
    assert_eq!(
        table[0],
        [
            "donor",
            "hydrogen",
            "acceptor",
            "donor_acceptor_distance_angstrom",
            "hydrogen_acceptor_distance_angstrom",
            "angle_degrees"
        ]
    );
    assert!(table.len() > 1, "a folded peptide has hydrogen bonds");
    let distance = column(&table, "donor_acceptor_distance_angstrom");
    let angle = column(&table, "angle_degrees");
    for row in &table[1..] {
        let d: f64 = row[distance].parse().expect("distance");
        let a: f64 = row[angle].parse().expect("angle");
        assert!(d <= 3.5 + 1e-4 && a >= 120.0 - 1e-4, "{row:?}");
    }
}

#[test]
fn a_stricter_angle_never_adds_hydrogen_bonds() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let count = |angle: &str| {
        invoke(
            directory.path(),
            &["interactions", "hbonds"],
            &[
                "--max-donor-acceptor-distance",
                "3.5",
                "--min-angle-degrees",
                angle,
            ],
        )
        .len()
    };
    assert!(count("150") <= count("120"));
}

#[test]
fn salt_bridges_are_anion_cation_pairs_within_the_distance() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = invoke(
        directory.path(),
        &["interactions", "salt-bridges"],
        &["--max-distance", "4.0"],
    );
    assert_eq!(table[0], ["anion", "cation", "distance_angstrom"]);
    let distance = column(&table, "distance_angstrom");
    for row in &table[1..] {
        let d: f64 = row[distance].parse().expect("distance");
        assert!(d <= 4.0 + 1e-4, "{row:?}");
    }
}

#[test]
fn stacking_names_its_arrangement() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = invoke(
        directory.path(),
        &["interactions", "pi-stacking"],
        &[
            "--max-centre-distance",
            "7.0",
            "--max-parallel-angle",
            "30",
            "--min-perpendicular-angle",
            "60",
            "--eigen-relative-tolerance",
            "1e-14",
            "--eigen-maximum-sweeps",
            "24",
        ],
    );
    assert_eq!(
        table[0],
        [
            "first_residue",
            "second_residue",
            "centre_distance_angstrom",
            "angle_degrees",
            "kind"
        ]
    );
    let kind = column(&table, "kind");
    for row in &table[1..] {
        assert!(
            ["parallel", "t-shaped"].contains(&row[kind].as_str()),
            "{row:?}"
        );
    }
}

#[test]
fn cation_pi_and_water_bridge_tables_have_their_documented_shape() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let cation = invoke(
        directory.path(),
        &["interactions", "cation-pi"],
        &[
            "--max-distance",
            "6.0",
            "--max-face-angle",
            "45",
            "--eigen-relative-tolerance",
            "1e-14",
            "--eigen-maximum-sweeps",
            "24",
        ],
    );
    assert_eq!(
        cation[0],
        ["cation_residue", "ring_residue", "distance_angstrom"]
    );
    let bridges = invoke(
        directory.path(),
        &["interactions", "water-bridges"],
        &[
            "--max-donor-acceptor-distance",
            "3.5",
            "--min-angle-degrees",
            "120",
        ],
    );
    assert_eq!(bridges[0], ["water", "first", "second"]);
    assert_eq!(
        bridges.len(),
        1,
        "a peptide without solvent has no water bridge"
    );
}

#[test]
fn every_threshold_is_required() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let output = common::run(
        directory.path(),
        &["interactions", "hbonds", peptide.to_str().expect("path")],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--max-donor-acceptor-distance"));
}
