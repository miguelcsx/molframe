//! End-to-end contact maps, native contacts and normal modes.

mod common;
use common::{haemoglobin_dimer, hydrogenated_peptide, rigidly_moved_peptide, rows, run};

#[test]
fn a_wider_cutoff_only_adds_residue_contacts_and_each_respects_its_cutoff() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path()).display().to_string();
    let at = |cutoff: &str| {
        rows(
            directory.path(),
            &[
                "contact-map",
                &peptide,
                "--cutoff",
                cutoff,
                "--min-separation",
                "2",
            ],
        )
    };
    let narrow = at("3.5");
    assert_eq!(
        narrow[0],
        ["first_residue", "second_residue", "min_distance_angstrom"]
    );
    let wide = at("5.0");
    assert!(narrow.len() > 1 && wide.len() > narrow.len());
    for row in &narrow[1..] {
        let first: u32 = row[0].parse().expect("residue");
        let second: u32 = row[1].parse().expect("residue");
        let distance: f64 = row[2].parse().expect("distance");
        assert!(second >= first + 2, "{row:?}");
        assert!(distance <= 3.5 + 1e-4, "{row:?}");
        assert!(
            wide.contains(row),
            "a contact at 3.5 is a contact at 5.0: {row:?}"
        );
    }
}

#[test]
fn a_rigid_motion_keeps_every_native_contact() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = hydrogenated_peptide(directory.path());
    let moved = rigidly_moved_peptide(directory.path());
    let table = rows(
        directory.path(),
        &[
            "native-contacts",
            native.to_str().expect("path"),
            moved.to_str().expect("path"),
            "--cutoff",
            "4.5",
            "--retention",
            "1.001",
        ],
    );
    assert_eq!(table[0], ["native", "kept", "fraction"]);
    assert_eq!(table[1][0], table[1][1], "every contact is kept: {table:?}");
    assert!(table[1][0].parse::<u32>().expect("count") > 100);
    let fraction: f64 = table[1][2].parse().expect("fraction");
    assert!((fraction - 1.0).abs() < 1e-12);
}

#[test]
fn moving_a_chain_away_breaks_its_native_contacts() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let apart = haemoglobin_dimer(directory.path(), "apart.pdb", 40.0);
    let table = rows(
        directory.path(),
        &[
            "native-contacts",
            native.to_str().expect("path"),
            apart.to_str().expect("path"),
            "--cutoff",
            "4.5",
            "--retention",
            "1.0",
        ],
    );
    let fraction: f64 = table[1][2].parse().expect("fraction");
    assert!(fraction > 0.3 && fraction < 0.99, "{table:?}");
}

fn modes(
    directory: &std::path::Path,
    model: &str,
    table: &str,
    extra: &[&str],
) -> Vec<Vec<String>> {
    let dimer = haemoglobin_dimer(directory, "native.pdb", 0.0)
        .display()
        .to_string();
    let mut arguments = vec![
        "nma",
        dimer.as_str(),
        "--network",
        model,
        "--sites",
        "chain A and name CA",
        "--contact-distance",
        "10",
        "--mode-count",
        "5",
        "--zero-mode-tolerance",
        "1e-8",
        "--memory-limit",
        "200000000",
        "--table",
        table,
    ];
    arguments.extend_from_slice(extra);
    rows(directory, &arguments)
}

#[test]
fn gnm_modes_are_positive_and_ordered_and_every_site_fluctuates() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = modes(
        directory.path(),
        "gnm",
        "modes",
        &["--reduction", "deterministic"],
    );
    assert_eq!(table[0], ["mode", "eigenvalue"]);
    assert_eq!(table.len(), 6);
    let values: Vec<f64> = table[1..]
        .iter()
        .map(|row| row[1].parse().expect("eigenvalue"))
        .collect();
    assert!(
        values[0] > 0.0 && values.windows(2).all(|pair| pair[0] <= pair[1]),
        "{values:?}"
    );
    let sites = modes(
        directory.path(),
        "gnm",
        "fluctuations",
        &["--reduction", "fast"],
    );
    assert_eq!(sites[0], ["site", "fluctuation"]);
    assert_eq!(
        sites.len(),
        142,
        "one row per alpha carbon of the 141-residue chain"
    );
    assert!(
        sites[1..]
            .iter()
            .all(|row| row[1].parse::<f64>().expect("fluctuation") > 0.0)
    );
}

#[test]
fn anm_has_six_zero_modes_and_gnm_requires_its_reduction() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let dimer = haemoglobin_dimer(directory.path(), "native.pdb", 0.0)
        .display()
        .to_string();
    let output = run(
        directory.path(),
        &[
            "nma",
            &dimer,
            "--network",
            "anm",
            "--sites",
            "chain A and name CA",
            "--contact-distance",
            "15",
            "--mode-count",
            "3",
            "--zero-mode-tolerance",
            "1e-6",
            "--memory-limit",
            "400000000",
            "--table",
            "modes",
            "--format",
            "tsv",
        ],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero modes: 6"));
    let refused = run(
        directory.path(),
        &[
            "nma",
            &dimer,
            "--network",
            "gnm",
            "--sites",
            "chain A and name CA",
            "--contact-distance",
            "10",
            "--mode-count",
            "3",
            "--zero-mode-tolerance",
            "1e-8",
            "--memory-limit",
            "400000000",
            "--table",
            "modes",
        ],
    );
    assert!(!refused.status.success());
    assert!(String::from_utf8_lossy(&refused.stderr).contains("--reduction"));
}
