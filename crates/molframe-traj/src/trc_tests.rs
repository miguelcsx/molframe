use super::*;

#[test]
fn shifted_and_rotated_boxes_are_retained() {
    let source = "TITLE\nexample\nEND\nTIMESTEP\n2 0.5\nEND\nPOSITIONRED\n0.1 0.2 0.3\nEND\nGENBOX\n1\n1 2 3\n90 90 90\n0.1 -0.2 0.3\n180 0 90\nEND\n";
    let parsed =
        parse_gromos11_trc(source).unwrap_or_else(|error| panic!("GROMOS parse failed: {error}"));
    assert_eq!(parsed.boundaries, [Some(GromosBoundary::Rectangular)]);
    assert_eq!(parsed.frames[0].positions, [[1.0, 2.0, 3.0]]);
    assert_eq!(
        parsed.frames[0].data.get("gromos11.origin_nm"),
        Some(&FrameValue::Floats(vec![0.1, -0.2, 0.3]))
    );
}

#[test]
fn topology_drift_is_rejected() {
    let source = "TITLE\nx\nEND\nPOSITIONRED\n0 0 0\nEND\nPOSITIONRED\n0 0 0\n1 1 1\nEND\n";
    assert!(matches!(
        parse_gromos11_trc(source),
        Err(GromosError::AtomCount {
            expected: 1,
            actual: 2
        })
    ));
}

#[test]
fn external_trc_fixture_when_configured() {
    let Some(path) = std::env::var_os("MOLFRAME_TRC_FIXTURE") else {
        return;
    };
    let source = std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let parsed =
        parse_gromos11_trc(&source).unwrap_or_else(|error| panic!("fixture parse failed: {error}"));
    assert!(!parsed.frames.is_empty());
    assert!(!parsed.frames[0].positions.is_empty());
}
