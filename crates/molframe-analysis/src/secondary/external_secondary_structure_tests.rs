use super::parse_dssp_output;

#[test]
fn annotated_mmcif_segments_are_typed() {
    let source = b"data_dssp\n\
loop_\n_struct_conf.id\n_struct_conf.conf_type_id\n_struct_conf.beg_label_asym_id\n\
_struct_conf.beg_label_seq_id\n_struct_conf.end_label_asym_id\n_struct_conf.end_label_seq_id\n\
1 HELX_RH_AL_P A 2 A 8\n2 STRN A 12 A 15\n";
    let segments =
        parse_dssp_output(source).unwrap_or_else(|error| panic!("DSSP output failed: {error}"));
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[0].kind.as_ref(), "HELX_RH_AL_P");
    assert_eq!(
        (segments[0].begin_sequence, segments[0].end_sequence),
        (2, 8)
    );
    assert_eq!(segments[1].kind.as_ref(), "STRN");
}

#[test]
fn malformed_output_is_reported() {
    assert!(parse_dssp_output(b"data_x\n_struct_conf.id 'unterminated\n").is_err());
}
