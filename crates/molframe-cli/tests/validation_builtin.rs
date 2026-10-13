//! End-to-end runs of the table-backed validation checks.

mod common;
use common::{bench_data, rows_with_status};

fn run(structure: &str, checks: &str, extra: &[&str]) -> (i32, Vec<Vec<String>>) {
    let directory = tempfile::tempdir().expect("scratch directory");
    let path = bench_data(structure).display().to_string();
    let mut arguments = vec!["validate", path.as_str(), "--checks", checks];
    arguments.extend_from_slice(extra);
    rows_with_status(directory.path(), &arguments)
}

fn value<'a>(rows: &'a [Vec<String>], check: &str, item: &str) -> &'a str {
    rows.iter()
        .find(|row| row[0] == check && row[1] == item)
        .map(|row| row[2].as_str())
        .expect("row present")
}

#[test]
fn ubiquitin_has_no_ramachandran_outliers_and_names_its_tables() {
    let (_, rows) = run("1ubq.cif", "protein", &[]);
    assert!(value(&rows, "ramachandran", "outlier").starts_with("0/"));
    assert!(value(&rows, "protein", "table:top8000-ramachandran").contains("CC-BY-4.0"));
}

#[test]
fn summary_rows_alone_do_not_fail_the_run() {
    let (status, _) = run("1crn.cif", "clashscore", &[]);
    assert_eq!(status, 0);
}

#[test]
fn a_missing_user_table_is_refused() {
    let (status, _) = run(
        "1crn.cif",
        "protein",
        &["--rama-table", "/nonexistent/rama.bin"],
    );
    assert_ne!(status, 0);
}
