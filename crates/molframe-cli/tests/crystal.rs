//! End-to-end crystallographic data commands on files built with known contents.

mod common;
use common::{bench_data, rows, rows_with_status};
use molframe::crystal::{
    DensityMap, ReflectionColumn, ReflectionColumnType, ReflectionDataset, ReflectionTable,
    ReflectionValue, write_mtz,
};
use molframe_core::structure::UnitCell;
use std::path::{Path, PathBuf};

fn mtz(directory: &Path) -> PathBuf {
    let integers = |items: &[i64]| {
        items
            .iter()
            .copied()
            .map(ReflectionValue::Integer)
            .collect()
    };
    let column = |label: &str, kind, values, dataset_id, code| ReflectionColumn {
        label: label.into(),
        column_type: kind,
        values,
        dataset_id,
        mtz_type: Some(code),
    };
    let table = ReflectionTable {
        title: "native".into(),
        cell: Some(UnitCell {
            lengths: [10.0, 11.0, 12.0],
            angles: [90.0, 91.0, 92.0],
        }),
        space_group_number: Some(1),
        space_group_name: Some("P 1".into()),
        columns: vec![
            column(
                "H",
                ReflectionColumnType::MillerIndex,
                integers(&[0, 1]),
                0,
                'H',
            ),
            column(
                "K",
                ReflectionColumnType::MillerIndex,
                integers(&[0, 0]),
                0,
                'H',
            ),
            column(
                "L",
                ReflectionColumnType::MillerIndex,
                integers(&[1, 0]),
                0,
                'H',
            ),
            column(
                "FP",
                ReflectionColumnType::Amplitude,
                vec![ReflectionValue::Real(10.5), ReflectionValue::Missing],
                1,
                'F',
            ),
        ],
        datasets: vec![ReflectionDataset {
            id: 1,
            project: "proj".into(),
            crystal: "xtal".into(),
            name: "native".into(),
            wavelength: Some(1.0),
            cell: None,
        }],
        history: vec!["created by molframe".into()],
        symmetry_operations: vec!["X,Y,Z".into()],
        sort_order: [1, 2, 3, 0, 0],
        resolution_range: None,
        missing_value: None,
        extra_header_records: Vec::new(),
    };
    let path = directory.join("native.mtz");
    std::fs::write(&path, write_mtz(&table).expect("encode")).expect("write mtz");
    path
}

fn map(directory: &Path, name: &str, values: Vec<f32>) -> PathBuf {
    let map = DensityMap {
        dimensions: [2, 2, 2],
        starts: [0; 3],
        sampling: [2, 2, 2],
        cell: UnitCell {
            lengths: [10.0; 3],
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values,
    };
    let path = directory.join(name);
    std::fs::write(&path, map.to_mrc_bytes().expect("encode")).expect("write map");
    path
}

fn property(table: &[Vec<String>], name: &str) -> String {
    table
        .iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("no {name} in {table:?}"))[1]
        .clone()
}

#[test]
fn an_mtz_file_reports_its_header_and_its_columns() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let file = mtz(directory.path()).display().to_string();
    let summary = rows(
        directory.path(),
        &["crystal", "mtz-info", &file, "--table", "summary"],
    );
    assert_eq!(summary[0], ["property", "value"]);
    assert_eq!(property(&summary, "title"), "native");
    assert_eq!(property(&summary, "reflections"), "2");
    assert_eq!(property(&summary, "space_group_name"), "P 1");
    assert_eq!(property(&summary, "cell_a"), "10");
    assert_eq!(property(&summary, "cell_gamma"), "92");
    let columns = rows(
        directory.path(),
        &["crystal", "mtz-info", &file, "--table", "columns"],
    );
    assert_eq!(columns[0], ["label", "type", "dataset"]);
    let labels: Vec<&str> = columns[1..].iter().map(|row| row[0].as_str()).collect();
    assert!(labels.contains(&"FP") && labels.contains(&"H"));
    let fp = columns
        .iter()
        .find(|row| row[0] == "FP")
        .expect("FP column");
    assert_eq!(fp[1], "amplitude");
}

#[test]
fn map_statistics_match_the_values_the_map_was_built_from() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let file = map(
        directory.path(),
        "ramp.mrc",
        (1..=8_u8).map(f32::from).collect(),
    )
    .display()
    .to_string();
    let table = rows(
        directory.path(),
        &[
            "crystal",
            "map-stats",
            &file,
            "--histogram-bins",
            "7",
            "--histogram-range",
            "1",
            "8",
        ],
    );
    assert_eq!(property(&table, "count"), "8");
    assert_eq!(property(&table, "minimum"), "1");
    assert_eq!(property(&table, "maximum"), "8");
    let mean: f64 = property(&table, "mean").parse().expect("mean");
    let sigma: f64 = property(&table, "sigma").parse().expect("sigma");
    assert!((mean - 4.5).abs() < 1e-9);
    assert!(
        (sigma - 5.25_f64.sqrt()).abs() < 1e-9,
        "population sigma of 1..=8"
    );
    let counted: usize = table
        .iter()
        .filter(|row| row[0].starts_with("histogram_bin_"))
        .map(|row| row[1].parse::<usize>().expect("count"))
        .sum();
    assert_eq!(counted, 8, "every voxel falls in some bin");
}

#[test]
fn a_linear_map_correlates_perfectly_and_its_reverse_perfectly_negatively() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let observed = map(
        directory.path(),
        "observed.mrc",
        (1..=8_u8).map(f32::from).collect(),
    );
    let scaled = map(
        directory.path(),
        "scaled.mrc",
        (1..=8_u8).map(|v| 2.0 * f32::from(v) + 1.0).collect(),
    );
    let reversed = map(
        directory.path(),
        "reversed.mrc",
        (1..=8_u8).rev().map(f32::from).collect(),
    );
    let coefficient = |other: &Path| -> f64 {
        let table = rows(
            directory.path(),
            &[
                "crystal",
                "map-correlation",
                observed.to_str().expect("path"),
                other.to_str().expect("path"),
            ],
        );
        property(&table, "coefficient")
            .parse()
            .expect("coefficient")
    };
    assert!((coefficient(&scaled) - 1.0).abs() < 1e-6);
    assert!((coefficient(&reversed) + 1.0).abs() < 1e-6);
}

#[test]
fn haemoglobin_has_symmetry_mates_within_the_stated_cutoff() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let file = bench_data("4hhb.pdb").display().to_string();
    let table = rows(
        directory.path(),
        &[
            "crystal",
            "mates",
            &file,
            "--cutoff",
            "4.0",
            "--candidate-limit",
            "4000000000",
        ],
    );
    assert_eq!(table[0][0], "source_atom");
    assert!(
        table.len() > 1,
        "crystal contacts exist in a packed crystal"
    );
    for row in &table[1..] {
        let distance: f64 = row[6].parse().expect("distance");
        assert!(distance <= 4.0 + 1e-6, "{row:?}");
        let lattice: [i32; 3] = [
            row[3].parse().expect("a"),
            row[4].parse().expect("b"),
            row[5].parse().expect("c"),
        ];
        assert!(
            lattice.iter().any(|&shift| shift != 0) || row[2] != "0",
            "a mate is a different copy, not the deposited one: {row:?}"
        );
    }
}

#[test]
fn a_structure_without_a_recorded_space_group_has_no_mates() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let peptide = common::hydrogenated_peptide(directory.path())
        .display()
        .to_string();
    let (status, table) = rows_with_status(
        directory.path(),
        &[
            "crystal",
            "mates",
            &peptide,
            "--cutoff",
            "4.0",
            "--candidate-limit",
            "1000",
        ],
    );
    assert_eq!(status, 7, "refused");
    assert!(
        table
            .iter()
            .all(|row| row.len() < 2 || row[0] != "source_atom")
    );
}
