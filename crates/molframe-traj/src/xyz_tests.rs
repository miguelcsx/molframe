use super::{parse_xyz, write_xyz};
use molframe_core::element::Element;

const TWO_FRAMES: &str = "2
frame one
C 0.0 0.0 0.0
O 1.2 0.0 0.0
2
frame two
C 0.0 0.0 1.0
O 1.2 0.0 1.0
";

#[test]
fn every_frame_is_read() {
    let Some(frames) = parse_xyz(TWO_FRAMES) else {
        panic!("valid XYZ");
    };
    assert_eq!(frames.len(), 2);
    assert_eq!(frames[0].comment, "frame one");
    assert_eq!(frames[0].atoms.len(), 2);
    assert_eq!(frames[0].atoms[0].element, Element::CARBON);
    assert!((frames[1].atoms[1].position[2] - 1.0).abs() < 1e-4);
}

#[test]
fn a_frame_converts_to_a_bare_coordinate_frame() {
    let Some(frames) = parse_xyz(TWO_FRAMES) else {
        panic!("valid");
    };
    let coordinate = frames[0].to_frame();
    assert_eq!(coordinate.positions.len(), 2);
    assert!(
        coordinate.positions[0]
            .iter()
            .all(|value| value.abs() < 1e-6)
    );
}

#[test]
fn a_truncated_frame_fails_to_parse() {
    assert!(parse_xyz("3\ncomment\nC 0 0 0\n").is_none());
}

#[test]
fn writer_preserves_elements_comments_and_coordinates() {
    let Some(frames) = parse_xyz(TWO_FRAMES) else {
        panic!("valid fixture");
    };
    assert_eq!(parse_xyz(&write_xyz(&frames)), Some(frames));
}
