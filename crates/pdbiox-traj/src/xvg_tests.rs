use super::{SeriesError, read_xvg};

const ENERGY: &str = r#"# This file was created by gmx energy
# Command line: gmx energy -f ener.edr
@    title "Gromacs Energies"
@    xaxis  label "Time (ps)"
@    yaxis  label "(kJ/mol)"
@ s0 legend "Potential"
@ s1 legend "Kinetic En."
    0.000000  -125000.5    21000.0
   10.000000  -125100.5    21050.0
   20.000000  -125200.5    21100.0
"#;

#[test]
fn legends_and_axis_labels_are_read_from_the_metadata_lines() {
    let Ok(series) = read_xvg(ENERGY.as_bytes()) else {
        panic!("the fixture table reads")
    };
    assert_eq!(series.title.as_deref(), Some("Gromacs Energies"));
    assert_eq!(series.abscissa_label.as_deref(), Some("Time (ps)"));
    assert_eq!(series.ordinate_label.as_deref(), Some("(kJ/mol)"));
    assert_eq!(series.legends.len(), 2);
    assert_eq!(series.legends[0].as_ref(), "Potential");
    assert_eq!(series.legends[1].as_ref(), "Kinetic En.");
}

#[test]
fn every_column_is_stored_contiguously_in_record_order() {
    let Ok(series) = read_xvg(ENERGY.as_bytes()) else {
        panic!("the fixture table reads")
    };
    assert_eq!(series.len(), 3);
    assert_eq!(series.columns.len(), 2);
    let Some(potential) = series.column_by_legend("Potential") else {
        panic!("the potential column is named")
    };
    assert!((potential[0] - -125_000.5).abs() < 1e-6);
    assert!((potential[2] - -125_200.5).abs() < 1e-6);
    assert!((series.abscissa[1] - 10.0).abs() < 1e-9);
}

#[test]
fn a_column_can_be_sampled_between_records() {
    let Ok(series) = read_xvg(ENERGY.as_bytes()) else {
        panic!("the fixture table reads")
    };
    let Some(midpoint) = series.sample(1, 15.0) else {
        panic!("the kinetic column samples")
    };
    assert!((midpoint - 21_075.0).abs() < 1e-6, "sampled {midpoint}");
}

#[test]
fn sampling_outside_the_recorded_range_clamps_to_the_end_records() {
    let Ok(series) = read_xvg(ENERGY.as_bytes()) else {
        panic!("the fixture table reads")
    };
    let (Some(before), Some(after)) = (series.sample(0, -5.0), series.sample(0, 500.0)) else {
        panic!("clamping still produces a value")
    };
    assert!((before - -125_000.5).abs() < 1e-6);
    assert!((after - -125_200.5).abs() < 1e-6);
}

#[test]
fn a_ragged_row_is_refused_rather_than_padded() {
    let ragged = ENERGY.replace(
        "   20.000000  -125200.5    21100.0",
        "   20.000000  -125200.5",
    );
    assert!(matches!(
        read_xvg(ragged.as_bytes()),
        Err(SeriesError::RaggedRow {
            row: 2,
            found: 2,
            expected: 3
        })
    ));
}

#[test]
fn a_non_numeric_token_reports_its_row_and_column() {
    let broken = ENERGY.replace(
        "   10.000000  -125100.5    21050.0",
        "   10.000000  wrong    21050.0",
    );
    assert!(matches!(
        read_xvg(broken.as_bytes()),
        Err(SeriesError::InvalidNumber { row: 1, column: 1 })
    ));
}

#[test]
fn dataset_separators_and_comments_do_not_become_records() {
    let separated = format!("{ENERGY}&\n# trailing comment\n   30.000000  -125300.5    21150.0\n");
    let Ok(series) = read_xvg(separated.as_bytes()) else {
        panic!("a separated table still reads")
    };
    assert_eq!(series.len(), 4);
}

#[test]
fn a_table_without_metadata_still_reads_its_numbers() {
    let bare = "0.0 1.0\n1.0 2.0\n";
    let Ok(series) = read_xvg(bare.as_bytes()) else {
        panic!("a bare table reads")
    };
    assert_eq!(series.len(), 2);
    assert_eq!(series.legends.len(), 1);
    assert!(series.legends[0].is_empty());
    assert!(series.title.is_none());
}
