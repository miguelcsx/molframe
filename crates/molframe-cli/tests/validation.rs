//! End-to-end validation checks with outcomes that are known in advance.

mod common;
use common::{
    amino_acid_ccd, atom_records, backbone_roles, hydrogenated_peptide, mirrored_peptide,
    rows_with_status,
};
use std::path::Path;

const CCD_VERSION: &str = "wwPDB-2026-10-03";

fn validate(directory: &Path, structure: &Path, check: &[&str]) -> (i32, Vec<Vec<String>>) {
    let ccd = amino_acid_ccd();
    let structure = structure.display().to_string();
    let ccd = ccd.display().to_string();
    let mut arguments = vec![
        "validate",
        structure.as_str(),
        "--ccd",
        ccd.as_str(),
        "--ccd-version",
        CCD_VERSION,
    ];
    arguments.extend_from_slice(check);
    rows_with_status(directory, &arguments)
}

fn findings(table: &[Vec<String>]) -> &[Vec<String>] {
    assert_eq!(table[0], ["check", "item", "value"]);
    &table[1..]
}

#[test]
fn a_mirror_image_inverts_every_chiral_centre_and_the_original_has_none() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let flags = ["--checks", "chirality", "--chirality-minimum-volume", "0.1"];
    let (status, original) = validate(
        directory.path(),
        &hydrogenated_peptide(directory.path()),
        &flags,
    );
    assert_eq!(status, 0, "{original:?}");
    assert!(findings(&original).is_empty());
    let (status, mirrored) = validate(
        directory.path(),
        &mirrored_peptide(directory.path()),
        &flags,
    );
    assert_eq!(status, 6, "findings exit with the consistency code");
    let rows = findings(&mirrored);
    assert!(!rows.is_empty());
    assert!(
        rows.iter()
            .all(|row| row[0] == "chirality" && row[2].contains("Inverted"))
    );
}

#[test]
fn every_peptide_bond_is_cis_under_a_threshold_of_one_hundred_eighty() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let roles = backbone_roles(directory.path());
    let roles = roles.display().to_string();
    let count = |threshold: &str| {
        let (_, table) = validate(
            directory.path(),
            &peptide,
            &[
                "--checks",
                "cis-peptide",
                "--role-profile",
                &roles,
                "--cis-threshold-degrees",
                threshold,
            ],
        );
        findings(&table).len()
    };
    let everything = count("180");
    assert!(
        everything > 10,
        "a twenty-residue peptide has many peptide bonds: {everything}"
    );
    assert!(count("30") < everything);
}

#[test]
fn a_covalently_consistent_peptide_has_no_overvalent_atom_and_no_ligand_bonds() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let (status, table) = validate(directory.path(), &peptide, &["--checks", "valence"]);
    assert_eq!((status, findings(&table).len()), (0, 0), "{table:?}");
    let (status, table) = validate(
        directory.path(),
        &peptide,
        &["--checks", "ligand", "--ligand-bond-tolerance", "0.3"],
    );
    assert_eq!((status, findings(&table).len()), (0, 0), "{table:?}");
}

#[test]
fn ramachandran_support_decides_who_is_an_outlier() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let roles = backbone_roles(directory.path()).display().to_string();
    let library = directory.path().join("rama.toml");
    std::fs::write(
        &library,
        "[library]\nid='uniform'\nversion='1'\n[[distribution]]\nname='all'\n\
         x_edges=[-180.0,0.0,180.0]\ny_edges=[-180.0,0.0,180.0]\nweights=[1.0,1.0,1.0,1.0]\n",
    )
    .expect("write library");
    let library = library.display().to_string();
    let count = |minimum: &str| {
        let (_, table) = validate(
            directory.path(),
            &peptide,
            &[
                "--checks",
                "ramachandran",
                "--role-profile",
                &roles,
                "--reference-library",
                &library,
                "--basin",
                "alpha-right=all",
                "--ramachandran-minimum-probability",
                minimum,
            ],
        );
        findings(&table).len()
    };
    // Every bin of the uniform grid holds a quarter of the mass.
    assert_eq!(count("0.2"), 0);
    assert!(count("0.3") > 10);
}

#[test]
fn tls_predicts_the_b_factor_of_every_atom_of_its_group() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let atoms = atom_records(&peptide);
    let zero = "[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]";
    let groups = |translation: &str| {
        let path = directory.path().join("tls.toml");
        std::fs::write(
            &path,
            format!(
                "[[group]]\nid='all'\nselection='all'\norigin=[0.0,0.0,0.0]\nt={translation}\nl={zero}\ns={zero}\n"
            ),
        )
        .expect("write groups");
        path.display().to_string()
    };
    let check = |translation: &str| {
        let path = groups(translation);
        let (_, table) = validate(
            directory.path(),
            &peptide,
            &[
                "--checks",
                "tls",
                "--tls-groups",
                &path,
                "--tls-max-deviation",
                "1.0",
                "--tls-symmetry-tolerance",
                "1e-9",
            ],
        );
        findings(&table).len()
    };
    // The file records B = 0, so a zero tensor predicts it exactly and a
    // non-zero translation tensor is wrong for every atom.
    assert_eq!(check(zero), 0);
    assert_eq!(check("[[0.1,0.0,0.0],[0.0,0.1,0.0],[0.0,0.0,0.1]]"), atoms);
}

#[test]
fn a_check_refuses_to_run_without_its_thresholds() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let output = common::run(
        directory.path(),
        &[
            "validate",
            peptide.to_str().expect("path"),
            "--checks",
            "ligand",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--ligand-bond-tolerance"));
}

#[test]
fn rotamer_support_flags_serine_torsions_in_the_sparse_bin() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let profile = directory.path().join("rotamers.toml");
    std::fs::write(
        &profile,
        "[profile]\nid='ser'\nversion='1'\n[[definition]]\ncomponent_id='SER'\nchi_index=1\n\
         atoms=['N','CA','CB','OG']\ndistribution='ser-chi1'\n",
    )
    .expect("write profile");
    let library = directory.path().join("chi.toml");
    std::fs::write(
        &library,
        "[library]\nid='two-bins'\nversion='1'\n[[distribution]]\nname='ser-chi1'\n\
         edges=[-180.0,0.0,180.0]\nweights=[1.0,3.0]\n",
    )
    .expect("write library");
    let (profile, library) = (profile.display().to_string(), library.display().to_string());
    let run = |minimum: &str| {
        let (_, table) = validate(
            directory.path(),
            &peptide,
            &[
                "--checks",
                "rotamer",
                "--rotamer-profile",
                &profile,
                "--reference-library",
                &library,
                "--rotamer-minimum-probability",
                minimum,
            ],
        );
        findings(&table).to_vec()
    };
    assert!(run("0.0").is_empty(), "nothing is below zero support");
    // Bins hold 25% and 75% of the mass, so a 50% floor flags only the sparse bin.
    for row in run("0.5") {
        assert!(row[2].contains("probability=0.25"), "{row:?}");
    }
    assert!(
        !run("0.9").is_empty(),
        "every serine torsion is below a 90% floor"
    );
    assert!(run("0.9").len() >= run("0.5").len());
}

#[test]
fn plane_restraints_flag_planes_beyond_the_stated_tolerance() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = hydrogenated_peptide(directory.path());
    let planes = directory.path().join("planes.toml");
    std::fs::write(
        &planes,
        "[[plane]]\nid='residue-two'\nselection='resid 2 and not element H'\n",
    )
    .expect("write planes");
    let planes = planes.display().to_string();
    let run = |tolerance: &str| {
        let (_, table) = validate(
            directory.path(),
            &peptide,
            &[
                "--checks",
                "plane-restraints",
                "--plane-restraints",
                &planes,
                "--planarity-tolerance",
                tolerance,
                "--plane-relative-tolerance",
                "1e-14",
                "--plane-maximum-sweeps",
                "24",
            ],
        );
        findings(&table).to_vec()
    };
    assert!(run("100").is_empty());
    let flagged = run("0");
    assert_eq!(flagged.len(), 1);
    assert_eq!(flagged[0][1], "residue-two");
}
