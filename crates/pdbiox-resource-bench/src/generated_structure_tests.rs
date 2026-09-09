use super::*;

#[test]
fn virtual_windows_match_the_exact_concatenated_source() {
    let context = ExecutionContext::default();
    let mut source = RepeatedSource::new(b"abc", b"XY", 9, 8, &context).expect("source");
    assert_eq!(source.len_hint(), Some(9));
    assert_eq!(source.window(0, 8).expect("first").bytes(), b"abcXYXYX");
    assert_eq!(source.window(2, 7).expect("crossing").bytes(), b"cXYXYXY");
    assert_eq!(source.window(9, 1).expect("end").bytes(), b"");
}

#[test]
fn virtual_length_and_windows_cross_the_u32_byte_boundary() {
    let context = ExecutionContext::default();
    let minimum = u64::from(u32::MAX) + 128;
    let mut source = RepeatedSource::new(b"h", b"row\n", minimum, 64, &context).expect("source");
    assert!(source.len_hint().is_some_and(|length| length >= minimum));
    let start = u64::from(u32::MAX) - 7;
    let window = source.window(start, 32).expect("large offset");
    assert_eq!(window.start(), start);
    assert_eq!(window.bytes().len(), 32);
}
