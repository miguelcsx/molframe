//! End-to-end explicit polymer role contracts using the rebuilt CLI and real crambin.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn cli(directory: &Path) -> Command {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let options = molframe::analysis::DsspOptions::default();
    let mut command = Command::new(env!("CARGO_BIN_EXE_molframe"));
    command
        .current_dir(directory)
        .env("XDG_CONFIG_HOME", directory)
        .arg("sse")
        .arg(root.join("../molframe-bench/data/1crn.cif"))
        .arg("--ccd")
        .arg(root.join("../molframe-chem/data/CCD-amino-acids.cif"))
        .args(["--ccd-version", "wwPDB-2026-10-03", "--format", "json"])
        .arg("--electrostatic-prefactor")
        .arg(options.electrostatic_prefactor.to_string())
        .arg(format!(
            "--hydrogen-bond-energy={}",
            options.hydrogen_bond_energy
        ))
        .arg("--amide-hydrogen-distance")
        .arg(options.amide_hydrogen_distance.to_string())
        .arg("--minimum-sequence-separation")
        .arg(options.minimum_sequence_separation.to_string())
        .arg("--helix-offset")
        .arg(options.helix_offset.to_string())
        .arg("--three-ten-offset")
        .arg(options.three_ten_offset.to_string())
        .arg("--pi-offset")
        .arg(options.pi_offset.to_string())
        .arg("--turn-offsets")
        .arg(options.turn_offsets.start().to_string())
        .arg(options.turn_offsets.end().to_string())
        .arg("--bend-angle-degrees")
        .arg(options.bend_angle_degrees.to_string());
    command
}

fn with_profile(directory: &Path, document: &serde_json::Value) -> Output {
    let profile = directory.join("roles.json");
    std::fs::write(&profile, document.to_string()).expect("write caller profile");
    cli(directory)
        .arg("--role-profile")
        .arg(profile)
        .output()
        .expect("run rebuilt CLI")
}

#[test]
fn real_crambin_with_complete_amino_acid_ccd_and_explicit_backbone_roles_is_evaluable() {
    let directory = tempfile::tempdir().expect("temporary CLI context");
    let profile = serde_json::json!({
        "profile_id": "wwpdb-amino-acid-backbone-2026-10-03",
        "rules": [
            {"atom_name": "N", "component_kind": 1, "role": 1},
            {"atom_name": "CA", "component_kind": 1, "role": 2},
            {"atom_name": "C", "component_kind": 1, "role": 4},
            {"atom_name": "O", "component_kind": 1, "role": 8}
        ]
    });
    let output = with_profile(directory.path(), &profile);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let document: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("JSON envelope");
    let rows = document["result"].as_array().expect("assignment rows");
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let (structure, _) = molframe::read_with_options(
        root.join("../molframe-bench/data/1crn.cif"),
        &molframe::ReadOptions::new(),
    )
    .expect("read real crambin");
    assert_eq!(rows.len(), structure.residue_count());
    let canonical = [
        "unknown",
        "coil",
        "alpha-helix",
        "3-10-helix",
        "pi-helix",
        "polyproline",
        "other-helix",
        "beta-bridge",
        "strand",
        "turn",
        "bend",
    ];
    for (index, row) in rows.iter().enumerate() {
        let residue = row["residue"].as_str().expect("residue index");
        let state = row["secondary_structure"]
            .as_str()
            .expect("canonical state");
        assert_eq!(residue.parse::<usize>().expect("residue index"), index);
        assert!(canonical.contains(&state), "noncanonical {state}");
        assert_ne!(
            state, "unknown",
            "complete heavy backbone at residue {index}"
        );
    }
    assert!(
        rows.iter()
            .any(|row| row["secondary_structure"] == "alpha-helix")
    );
    assert!(
        rows.iter()
            .any(|row| row["secondary_structure"] == "strand")
    );
}

#[test]
fn missing_and_invalid_profiles_fail_without_emitting_an_assignment() {
    let directory = tempfile::tempdir().expect("temporary CLI context");
    let output = cli(directory.path())
        .output()
        .expect("run CLI without profile");
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("--role-profile"));
    assert!(output.stdout.is_empty());
    let output = cli(directory.path())
        .arg("--role-profile")
        .arg(directory.path().join("absent.toml"))
        .output()
        .expect("run CLI with missing file");
    assert_eq!(output.status.code(), Some(8));
    assert!(output.stdout.is_empty());
    for profile in [
        serde_json::json!({"profile_id": "v1", "rules": []}),
        serde_json::json!({"profile_id": "", "rules": [{"atom_name": "N", "component_kind": 1, "role": 1}]}),
        serde_json::json!({"profile_id": "v1", "rules": [{"atom_name": "N", "role": 1}]}),
        serde_json::json!({"profile_id": "v1", "rules": [{"atom_name": "", "component_kind": 1, "role": 1}]}),
        serde_json::json!({"profile_id": "v1", "rules": [{"atom_name": "N", "component_kind": 1, "role": 0}]}),
        serde_json::json!({"profile_id": "v1", "rules": [{"atom_name": "N", "component_kind": 1, "role": 262_144}]}),
        serde_json::json!({"profile_id": "v1", "rules": [{"atom_name": "N", "component_kind": 8, "role": 1}]}),
    ] {
        let output = with_profile(directory.path(), &profile);
        assert_eq!(output.status.code(), Some(8), "accepted {profile}");
        assert!(output.stdout.is_empty(), "assignment for invalid {profile}");
        assert!(!output.stderr.is_empty());
    }
}
