//! End-to-end mapped `DockQ` on real haemoglobin chains with the complete
//! amino-acid dictionary: renaming the chains must not change the score.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

const SHIFT_X: f64 = 3.0;
const SHIFT_Y: f64 = 1.5;

/// The alpha (A) and beta (B) chains of 4HHB, optionally renamed, with the beta
/// chain displaced when `displaced` so the result is a non-native pose.
fn dimer(rename: Option<(char, char)>, displaced: bool) -> String {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let source = std::fs::read_to_string(root.join("../molframe-bench/data/4hhb.pdb"))
        .expect("bundled 4HHB");
    let mut out = String::new();
    for line in source.lines().filter(|line| line.starts_with("ATOM")) {
        let chain = line.as_bytes()[21] as char;
        if chain != 'A' && chain != 'B' {
            continue;
        }
        let mut x: f64 = line[30..38].trim().parse().expect("x column");
        let mut y: f64 = line[38..46].trim().parse().expect("y column");
        if displaced && chain == 'B' {
            x += SHIFT_X;
            y += SHIFT_Y;
        }
        let label = match rename {
            Some((alpha, _)) if chain == 'A' => alpha,
            Some((_, beta)) => beta,
            None => chain,
        };
        let _ = writeln!(
            out,
            "{}{label}{}{x:8.3}{y:8.3}{}",
            &line[..21],
            &line[22..30],
            &line[46..]
        );
    }
    out.push_str("END\n");
    out
}

fn dockq(directory: &Path, model: &str, native: &str, mapped: bool) -> Vec<(String, f64)> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut command = Command::new(env!("CARGO_BIN_EXE_molframe"));
    command
        .current_dir(directory)
        .env("XDG_CONFIG_HOME", directory)
        .args(["compare", model, native, "--metrics", "dockq"])
        .args(["--receptor", "A", "--ligand", "B"])
        .args(["--contact-distance", "5", "--ligand-scale", "8.5"])
        .args(["--interface-scale", "1.5", "--format", "tsv"]);
    if mapped {
        command
            .arg("--map-chains")
            .arg("--ccd")
            .arg(root.join("../molframe-chem/data/CCD-amino-acids.cif"))
            .args(["--ccd-version", "wwPDB-2026-10-03"])
            .args(["--min-identity", "0.9", "--match-score", "2"])
            .args(["--mismatch-score=-1", "--gap-open=-5", "--gap-extend=-1"])
            .args(["--automorphism-limit", "64"]);
    }
    let output = command.output().expect("run the CLI");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .expect("utf-8 output")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.starts_with("metric"))
        .filter_map(|line| {
            let (name, value) = line.split_once('\t')?;
            Some((name.to_owned(), value.parse().ok()?))
        })
        .collect()
}

#[test]
fn renaming_the_chains_of_a_real_dimer_leaves_the_mapped_dockq_unchanged() {
    let directory = tempfile::tempdir().expect("temporary CLI context");
    let native = directory.path().join("native.pdb");
    let model = directory.path().join("model.pdb");
    let renamed = directory.path().join("renamed.pdb");
    std::fs::write(&native, dimer(None, false)).expect("write native");
    std::fs::write(&model, dimer(None, true)).expect("write model");
    std::fs::write(&renamed, dimer(Some(('Y', 'X')), true)).expect("write renamed model");
    let native = native.to_str().expect("utf-8 path");
    let plain = dockq(
        directory.path(),
        model.to_str().expect("path"),
        native,
        false,
    );
    let through_mapping = dockq(
        directory.path(),
        renamed.to_str().expect("path"),
        native,
        true,
    );
    assert_eq!(plain.len(), 4);
    assert_eq!(plain, through_mapping);
    let score = plain[0].1;
    assert!(
        score > 0.0 && score < 1.0,
        "a displaced pose scores between 0 and 1: {score}"
    );
}
