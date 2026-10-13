//! `molframe data` over a manifest, with no network: listing, verifying and the
//! errors a bad manifest or an absent entry give.

use std::path::Path;
use std::process::{Command, Output};

const MANIFEST: &str = r#"
[[entry]]
name = "lab-table"
version = "3"
url = "https://example.org/lab.tsv"
licence = "CC0-1.0"
sha256 = "1111111111111111111111111111111111111111111111111111111111111111"
extension = "lab"
converter = "tsv"
"#;

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_molframe"))
        .arg("data")
        .args(arguments)
        .output()
        .expect("the binary runs")
}

fn setup(directory: &Path) -> (String, String) {
    let manifest = directory.join("manifest.toml");
    std::fs::write(&manifest, MANIFEST).expect("manifest");
    let cache = directory.join("cache");
    std::fs::create_dir(&cache).expect("cache");
    (manifest.display().to_string(), cache.display().to_string())
}

#[test]
fn list_and_verify_report_a_missing_then_a_current_entry() {
    let directory = tempfile::tempdir().expect("scratch");
    let (manifest, cache) = setup(directory.path());
    let listed = run(&["list", "--manifest", &manifest, "--dir", &cache]);
    assert!(listed.status.success());
    let text = String::from_utf8_lossy(&listed.stdout);
    assert!(text.contains(r#""name":"lab-table""#) && text.contains(r#""status":"missing""#));
    assert!(
        !run(&["verify", "--manifest", &manifest, "--dir", &cache])
            .status
            .success()
    );
    let stored = format!(
        "# SPDX-License-Identifier: CC0-1.0\nformat: 1\nname: lab-table\nversion: 3\n\
         source: lab\nlicence: CC0-1.0\n# source-sha256: {}\n\nC\t1.7\n",
        "1".repeat(64)
    );
    std::fs::write(Path::new(&cache).join("lab-table.lab"), stored).expect("table");
    let verified = run(&[
        "verify",
        "lab-table",
        "--manifest",
        &manifest,
        "--dir",
        &cache,
    ]);
    assert!(verified.status.success());
    assert!(String::from_utf8_lossy(&verified.stdout).contains("current"));
}

#[test]
fn an_unknown_entry_and_an_invalid_manifest_are_refused() {
    let directory = tempfile::tempdir().expect("scratch");
    let (manifest, cache) = setup(directory.path());
    let unknown = run(&[
        "verify",
        "nothing",
        "--manifest",
        &manifest,
        "--dir",
        &cache,
    ]);
    assert_eq!(unknown.status.code(), Some(2));
    std::fs::write(&manifest, MANIFEST.replace("https://", "http://")).expect("manifest");
    let invalid = run(&["list", "--manifest", &manifest, "--dir", &cache]);
    assert_eq!(invalid.status.code(), Some(5));
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("https"));
}
