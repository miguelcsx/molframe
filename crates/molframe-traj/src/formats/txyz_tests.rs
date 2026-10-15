use super::{parse_txyz_records, write_txyz};

const FRAME: &str = "2  water\n\
1  O   0.0 0.0 0.0  1  2\n\
2  H   1.0 0.0 0.0  2  1\n";

#[test]
fn coordinates_are_read_past_the_index_and_name() {
    let frames = parse_txyz_records(FRAME).expect("valid TXYZ");
    assert_eq!(frames.len(), 1);
    assert_eq!(frames[0].atoms.len(), 2);
    assert!((frames[0].atoms[1].position[0] - 1.0).abs() < 1e-6);
}

#[test]
fn an_arc_of_two_frames_is_read() {
    let frames = parse_txyz_records(&format!("{FRAME}{FRAME}")).expect("valid ARC");
    assert_eq!(frames.len(), 2);
}

#[test]
fn a_short_frame_fails_to_parse() {
    assert!(parse_txyz_records("3  x\n1 O 0 0 0 1\n").is_err());
}

#[test]
fn written_arc_round_trips_every_frame() {
    let frames = parse_txyz_records(&format!("{FRAME}{FRAME}")).expect("complete ARC");
    let encoded = write_txyz(&frames).expect("write ARC");
    assert_eq!(parse_txyz_records(&encoded).expect("reparse ARC"), frames);
}
