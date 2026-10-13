//! End-to-end `sasa --radii-file`: a caller table replaces the set, reproduces the
//! default when it only names a base, and the applied policy is printed.

mod common;
use common::{bench_data, run};

fn sasa(directory: &std::path::Path, extra: &[&str]) -> std::process::Output {
    let input = bench_data("1crn.cif").display().to_string();
    let mut arguments = vec![
        "sasa", &input, "--probe", "1.4", "--points", "100", "--radii", "bondi",
    ];
    arguments.extend_from_slice(extra);
    run(directory, &arguments)
}

#[test]
fn a_base_only_table_matches_the_default_and_the_policy_is_reported() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = directory.path().join("radii.txt");
    std::fs::write(
        &table,
        "name: lab\nversion: 3\nsource: bench notes\nlicence: CC0\nbase: bondi\n",
    )
    .expect("write table");
    let default = sasa(directory.path(), &[]);
    let custom = sasa(
        directory.path(),
        &["--radii-file", &table.display().to_string()],
    );
    assert!(default.status.success() && custom.status.success());
    let areas = |output: &std::process::Output| {
        let document: serde_json::Value =
            serde_json::from_slice(&output.stdout).expect("JSON output");
        document["result"].clone()
    };
    assert_eq!(areas(&default), areas(&custom));
    let header = String::from_utf8_lossy(&custom.stderr);
    assert!(header.contains("# radius policy:"), "{header}");
    assert!(
        header.contains("lab 3 (source: bench notes; licence: CC0)"),
        "{header}"
    );
    assert!(String::from_utf8_lossy(&default.stderr).contains("bondi"));
}

#[test]
fn a_malformed_table_is_a_usage_error_naming_the_line() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let table = directory.path().join("radii.txt");
    std::fs::write(&table, "name: lab\nC -1\n").expect("write table");
    let output = sasa(
        directory.path(),
        &["--radii-file", &table.display().to_string()],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("line"));
}
