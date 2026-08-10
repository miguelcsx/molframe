use super::*;

const LINE: &str = "ATOM      1  N   GLU L   1      27.340  24.430   2.614  1.00  0.00           N";

#[test]
fn fields_are_one_based_inclusive_and_trimmed() {
    assert_eq!(record(LINE), "ATOM");
    assert_eq!(text(LINE, 13, 16), "N");
    assert_eq!(text(LINE, 18, 20), "GLU");
    assert_eq!(text(LINE, 22, 22), "L");
    assert_eq!(integer(LINE, 7, 11), Some(1));
    assert_eq!(integer(LINE, 23, 26), Some(1));
    assert_eq!(real(LINE, 31, 38), Some(27.340));
    assert_eq!(text(LINE, 77, 78), "N");
}

#[test]
fn an_untrimmed_field_keeps_the_leading_space_that_names_the_element() {
    assert_eq!(raw(LINE, 13, 16), " N  ");
    assert_eq!(raw("HETATM    1 ZN    ZN A 100", 13, 16), "ZN  ");
}

#[test]
fn a_line_that_stops_short_reads_as_though_it_were_blank_padded() {
    let short = "ATOM      1  N   GLU L   1";
    assert_eq!(text(short, 31, 38), "");
    assert_eq!(real(short, 31, 38), None);
    assert_eq!(text(short, 13, 16), "N");
    assert_eq!(raw(short, 77, 78), "");
}

#[test]
fn an_unparseable_field_yields_nothing_rather_than_a_wrong_number() {
    let broken = "ATOM    abc  N   GLU L   1";
    assert_eq!(integer(broken, 7, 11), None);
    assert_eq!(real(broken, 7, 11), None);
}

#[test]
fn an_empty_line_yields_empty_fields_at_every_position() {
    assert_eq!(record(""), "");
    assert_eq!(text("", 1, 80), "");
    assert_eq!(integer("", 7, 11), None);
}
