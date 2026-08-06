use super::*;

#[test]
fn advancing_over_a_newline_starts_the_next_line_at_column_one() {
    let at = Position::START.advance(b'a').advance(b'\n');
    assert_eq!(at, Position::new(2, 2, 1));
}

#[test]
fn advancing_over_an_ordinary_byte_moves_the_column_only() {
    let at = Position::START.advance(b'a').advance(b'b');
    assert_eq!(at, Position::new(2, 1, 3));
}

#[test]
fn a_span_yields_exactly_the_bytes_it_covers() {
    let source = b"data_x\nloop_\n";
    let span = ByteSpan::new(Position::new(7, 2, 1), 12);
    assert_eq!(span.slice(source), Some(&b"loop_"[..]));
    assert_eq!(span.len(), 5);
}

#[test]
fn a_span_past_the_end_of_its_buffer_refuses_rather_than_clamping() {
    let source = b"short";
    let span = ByteSpan::new(Position::new(0, 1, 1), 99);
    assert_eq!(span.slice(source), None);
}

#[test]
fn an_empty_span_marks_a_point_without_covering_a_byte() {
    let span = ByteSpan::empty(Position::new(3, 1, 4));
    assert!(span.is_empty());
    assert_eq!(span.slice(b"abcdef"), Some(&b""[..]));
}
