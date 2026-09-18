use super::*;

#[test]
fn advancing_over_a_newline_starts_the_next_line_at_column_one() {
    let at = Position::START
        .advance(b'a')
        .and_then(|position| position.advance(b'\n'));
    assert_eq!(at, Some(Position::new(2, 2, 1)));
}

#[test]
fn advancing_over_an_ordinary_byte_moves_the_column_only() {
    let at = Position::START
        .advance(b'a')
        .and_then(|position| position.advance(b'b'));
    assert_eq!(at, Some(Position::new(2, 1, 3)));
}

#[test]
fn a_span_yields_exactly_the_bytes_it_covers() {
    let source = b"data_x\nloop_\n";
    let span = ByteSpan::new(Position::new(7, 2, 1), 12);
    assert_eq!(span.slice(source), Some(&b"loop_"[..]));
    assert_eq!(span.len(), Some(5));
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

#[test]
fn a_position_addresses_offsets_past_four_gibibytes() {
    let beyond = u64::from(u32::MAX) + 1;
    let at = Position::new(beyond, 1, 1);
    assert_eq!(at.byte_offset, beyond);

    let Some(next) = at.advance(b'x') else {
        panic!("advancing past four gibibytes must not overflow")
    };
    assert_eq!(next.byte_offset, beyond + 1);
    assert_eq!(next.column, 2);
}

#[test]
fn a_span_covers_a_range_past_four_gibibytes() {
    let start = Position::new(u64::from(u32::MAX) + 1, 7, 1);
    let span = ByteSpan::new(start, start.byte_offset + 12);
    assert_eq!(span.len(), Some(12));
    assert!(!span.is_empty());
}

#[test]
fn findings_order_by_offset_across_the_four_gibibyte_boundary() {
    let early = ByteSpan::empty(Position::new(u64::from(u32::MAX) - 1, 1, 1));
    let late = ByteSpan::empty(Position::new(u64::from(u32::MAX) + 1, 1, 1));
    assert!(early.start.byte_offset < late.start.byte_offset);
}
