//! End-to-end comparison metrics with outcomes known in advance.

mod common;
use common::{
    amino_acid_ccd, haemoglobin_dimer, hydrogenated_peptide, rigidly_moved_peptide,
    ring_flipped_peptide, rows,
};
use std::path::Path;

/// Runs `compare MODEL REFERENCE` and returns each printed metric as (name, value).
fn compare(directory: &Path, model: &Path, reference: &Path, extra: &[&str]) -> Vec<(String, f64)> {
    let (model, reference) = (model.display().to_string(), reference.display().to_string());
    let mut arguments = vec!["compare", model.as_str(), reference.as_str()];
    arguments.extend_from_slice(extra);
    rows(directory, &arguments)
        .into_iter()
        .skip(1)
        .map(|row| (row[0].clone(), row[1].parse().expect("numeric value")))
        .collect()
}

fn close(left: f64, right: f64) -> bool {
    (left - right).abs() < 1e-9
}

fn value(table: &[(String, f64)], name: &str) -> f64 {
    table
        .iter()
        .find(|(candidate, _)| candidate == name)
        .unwrap_or_else(|| panic!("no metric {name} in {table:?}"))
        .1
}

#[test]
fn a_rigid_motion_disappears_under_a_fit_and_shows_without_one() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = hydrogenated_peptide(directory.path());
    let moved = rigidly_moved_peptide(directory.path());
    let region = ["--metrics", "interface-rmsd", "--region", "resid 1:8"];
    let fitted = compare(
        directory.path(),
        &moved,
        &native,
        &[&region[..], &["--alignment", "fit"]].concat(),
    );
    assert!(value(&fitted, "interface_rmsd") < 1e-2, "{fitted:?}");
    let unfitted = compare(
        directory.path(),
        &moved,
        &native,
        &[&region[..], &["--alignment", "none"]].concat(),
    );
    assert!(value(&unfitted, "interface_rmsd") > 1.0, "{unfitted:?}");
    let pocket = compare(
        directory.path(),
        &moved,
        &native,
        &[
            "--metrics",
            "pocket-rmsd",
            "--region",
            "resid 1:8",
            "--alignment",
            "fit",
        ],
    );
    assert!(value(&pocket, "pocket_rmsd") < 1e-2);
}

#[test]
fn a_ring_flip_costs_nothing_under_chemical_symmetry() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = hydrogenated_peptide(directory.path());
    let (flipped, tyrosine) = ring_flipped_peptide(directory.path());
    let ccd = amino_acid_ccd().display().to_string();
    let table = compare(
        directory.path(),
        &flipped,
        &native,
        &[
            "--metrics",
            "ligand-rmsd",
            "--ccd",
            &ccd,
            "--ccd-version",
            "wwPDB-2026-10-03",
            "--ligand-selection",
            &format!("resname TYR and resid {tyrosine}"),
            "--ligand-component",
            "TYR",
            "--ligand-automorphism-limit",
            "64",
        ],
    );
    assert!(value(&table, "ligand_rmsd") < 1e-5, "{table:?}");
    // Without the symmetry the swapped halves of the ring would not coincide.
    let (flipped, tyrosine) = ring_flipped_peptide(directory.path());
    let plain = compare(
        directory.path(),
        &flipped,
        &native,
        &[
            "--metrics",
            "pocket-rmsd",
            "--region",
            &format!("resname TYR and resid {tyrosine}"),
            "--alignment",
            "none",
        ],
    );
    assert!(value(&plain, "pocket_rmsd") > 0.5, "{plain:?}");
}

#[test]
fn identical_dimers_agree_perfectly_on_every_contact_score() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let same = haemoglobin_dimer(directory.path(), "same.pdb", 0.0);
    let table = compare(
        directory.path(),
        &same,
        &native,
        &[
            "--metrics",
            "qs,contact-similarity",
            "--receptor",
            "A",
            "--ligand",
            "B",
            "--contact-distance",
            "5",
            "--min-separation",
            "1",
        ],
    );
    assert!((value(&table, "qs") - 1.0).abs() < 1e-12);
    assert!((value(&table, "contact_similarity") - 1.0).abs() < 1e-12);
    assert!(close(
        value(&table, "contact_similarity_shared"),
        value(&table, "contact_similarity_union")
    ));
}

#[test]
fn separating_the_chains_loses_interface_contacts() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let apart = haemoglobin_dimer(directory.path(), "apart.pdb", 40.0);
    let table = compare(
        directory.path(),
        &apart,
        &native,
        &[
            "--metrics",
            "qs",
            "--receptor",
            "A",
            "--ligand",
            "B",
            "--contact-distance",
            "5",
        ],
    );
    assert!(value(&table, "qs") < 1e-12, "{table:?}");
}

#[test]
fn identical_structures_have_zero_lost_contact_area() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = hydrogenated_peptide(directory.path());
    let table = compare(
        directory.path(),
        &native,
        &native,
        &[
            "--metrics",
            "cad",
            "--cad-probe",
            "1.4",
            "--cad-density",
            "2.0",
            "--radii",
            "bondi",
        ],
    );
    assert!(value(&table, "cad_lost_area").abs() < 1e-9);
    assert!(value(&table, "cad_reference_area") > 0.0);
    assert!((value(&table, "cad") - 1.0).abs() < 1e-12);
}

#[test]
fn ce_aligns_a_chain_with_itself_atom_for_atom() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = haemoglobin_dimer(directory.path(), "native.pdb", 0.0);
    let table = compare(
        directory.path(),
        &native,
        &native,
        &[
            "--metrics",
            "ce",
            "--guide",
            "chain A and name CA",
            "--ce-window",
            "8",
            "--ce-max-gap",
            "30",
            "--ce-max-paths",
            "20",
            "--ce-fragment-threshold=-3.0",
            "--ce-path-threshold=-4.0",
            "--ce-memory-limit",
            "100000000",
        ],
    );
    assert!(value(&table, "ce_rmsd") < 1e-3, "{table:?}");
    // CE aligns whole fragments of eight guide atoms: the 141 alpha-chain
    // residues hold seventeen of them, and the last five cannot form one.
    assert!(close(value(&table, "ce_fragments"), 17.0));
    assert!(close(value(&table, "ce_aligned_atoms"), 136.0));
}

#[test]
fn a_metric_refuses_to_run_without_its_parameters() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let native = hydrogenated_peptide(directory.path());
    let output = common::run(
        directory.path(),
        &[
            "compare",
            native.to_str().expect("path"),
            native.to_str().expect("path"),
            "--metrics",
            "cad",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--cad-probe"));
}
