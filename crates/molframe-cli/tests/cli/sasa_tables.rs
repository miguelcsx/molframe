//! End-to-end `sasa --radii <PATH>`: a caller table replaces the profile's, a copy
//! of the default reproduces it, and the applied table is printed.

use crate::common::{bench_data, run};

fn sasa(directory: &std::path::Path, radii: &str) -> std::process::Output {
    let input = bench_data("1crn.cif").display().to_string();
    run(
        directory,
        &[
            "sasa", &input, "--probe", "1.4", "--points", "100", "--radii", radii,
        ],
    )
}

#[test]
fn a_copy_of_the_default_matches_it_and_the_applied_table_is_reported() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = directory.path().join("lab.radii");
    let bondi = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/bondi.radii"),
    )
    .expect("bondi table");
    std::fs::write(
        &table,
        bondi
            .replace("name: bondi", "name: lab")
            .replace("version: 1964", "version: 3")
            .replace("licence: MIT", "licence: CC0"),
    )
    .expect("write table");
    let default = sasa(directory.path(), "bondi");
    let custom = sasa(directory.path(), &table.display().to_string());
    assert!(default.status.success() && custom.status.success());
    let areas = |output: &std::process::Output| {
        let document: serde_json::Value =
            serde_json::from_slice(&output.stdout).expect("JSON output");
        document["result"].clone()
    };
    assert_eq!(areas(&default), areas(&custom));
    let header = String::from_utf8_lossy(&custom.stderr);
    assert!(header.contains("# radius policy:"), "{header}");
    assert!(header.contains("lab"), "{header}");
    assert!(String::from_utf8_lossy(&default.stderr).contains("bondi"));
}

#[test]
fn a_malformed_table_is_a_usage_error_naming_the_line() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = directory.path().join("radii.radii");
    std::fs::write(
        &table,
        "format: 1\nname: lab\nversion: 1\nsource: s\nlicence: CC0\nC -1\n",
    )
    .expect("write table");
    let output = sasa(directory.path(), &table.display().to_string());
    assert!(!output.status.success());
    let text = String::from_utf8_lossy(&output.stderr);
    assert!(text.contains("E2201") && text.contains("line"), "{text}");
}
