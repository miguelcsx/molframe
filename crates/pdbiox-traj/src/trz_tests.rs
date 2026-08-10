use super::*;
use crate::write_trz;

#[test]
fn external_trz_fixture_when_configured() {
    let Some(path) = std::env::var_os("PDBIOX_TRZ_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let parsed = parse_trz(&bytes).unwrap_or_else(|error| panic!("fixture parse failed: {error}"));
    assert!(!parsed.frames.is_empty());
    assert_eq!(parsed.frames[0].positions.len(), 8184);
}

#[test]
fn parsed_frame_round_trips_without_losing_metadata() {
    let Some(path) = std::env::var_os("PDBIOX_TRZ_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let parsed = parse_trz(&bytes).unwrap_or_else(|error| panic!("fixture parse failed: {error}"));
    let rewritten = write_trz(&parsed).unwrap_or_else(|error| panic!("TRZ write failed: {error}"));
    let reparsed =
        parse_trz(&rewritten).unwrap_or_else(|error| panic!("TRZ reparse failed: {error}"));
    assert_eq!(reparsed.title, parsed.title);
    assert_eq!(reparsed.frames.len(), parsed.frames.len());
    assert!(
        reparsed.frames[0]
            .positions
            .iter()
            .flatten()
            .zip(parsed.frames[0].positions.iter().flatten())
            .all(|(left, right)| (*left - *right).abs() < 2.0e-5)
    );
    assert_eq!(reparsed.frames[0].data, parsed.frames[0].data);
}
