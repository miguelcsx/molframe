use super::{parse_charmm_record, write_charmm_card};

#[test]
fn standard_writer_round_trips_fixed_width_negative_coordinates() {
    let source = "* title\n*\n    1\n    1    1 MOL  CA  -123.45600   2.50000  99.00000 SYS  1      0.00000\n";
    let record = parse_charmm_record(source).expect("CARD record");
    let encoded = write_charmm_card(&record).expect("CARD write");
    assert_eq!(parse_charmm_record(&encoded).expect("CARD reparse"), record);
}

#[test]
fn extended_card_uses_expanded_coordinate_columns() {
    let atom = format!(
        "{:>10}{:>10}  {:<8}  {:<8}{:>20.10}{:>20.10}{:>20.10}  {:<8}  {:<8}{:>20.10}",
        1, 1, "MOL", "X", 1.25, -2.5, 3.75, "SYS", "1", 0.0,
    );
    let source = format!("* extended\n*\n         1 EXT\n{atom}\n");
    let frame = parse_charmm_record(&source)
        .unwrap_or_else(|error| panic!("extended CARD failed: {error}"))
        .to_frame();
    assert!((frame.positions[0][1] + 2.5).abs() < 1.0e-6);
}

#[test]
fn free_field_card_is_supported_explicitly() {
    let source = "* free\n1 FREE\n1 1 MOL X 1.0 2.0 3.0 SYS 1 0.0\n";
    let frame = parse_charmm_record(source)
        .unwrap_or_else(|error| panic!("free CARD failed: {error}"))
        .to_frame();
    assert!((frame.positions[0][2] - 3.0).abs() < f32::EPSILON);
}
