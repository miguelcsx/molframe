use super::*;

#[test]
fn both_endian_variants_round_trip() {
    let frame = Timestep {
        positions: vec![[1.25, -2.5, 3.75], [4.0, 5.0, -6.0]],
        ..Timestep::default()
    };
    for endian in [NamdEndian::Little, NamdEndian::Big] {
        let bytes = write_namd_binary(&frame, endian)
            .unwrap_or_else(|error| panic!("NAMD write failed: {error}"));
        let parsed =
            parse_namd_binary(&bytes).unwrap_or_else(|error| panic!("NAMD parse failed: {error}"));
        assert_eq!(parsed.endian, endian);
        assert_eq!(parsed.frame.positions, frame.positions);
    }
}

#[test]
fn exact_length_is_required() {
    assert_eq!(
        parse_namd_binary(&1_i32.to_le_bytes()),
        Err(NamdError::InvalidData)
    );
}

#[test]
fn external_namd_fixture_when_configured() {
    let Some(path) = std::env::var_os("MOLFRAME_NAMD_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let parsed =
        parse_namd_binary(&bytes).unwrap_or_else(|error| panic!("fixture parse failed: {error}"));
    assert_eq!(parsed.frame.positions.len(), 3341);
    assert_eq!(parsed.endian, NamdEndian::Little);
}
