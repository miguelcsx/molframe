use super::command::value_text;
use molframe::formats::cif::CifValue;

#[test]
fn cif_sentinels_remain_distinct_in_show_output() {
    assert_eq!(value_text(Some(&CifValue::Inapplicable)), ".");
    assert_eq!(value_text(Some(&CifValue::Unknown)), "?");
    assert_eq!(value_text(None), "?");
}
