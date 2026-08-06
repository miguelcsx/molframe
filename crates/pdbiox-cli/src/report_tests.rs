use super::*;

#[test]
fn an_object_renders_its_fields_in_the_order_they_were_added() {
    let mut object = Json::new();
    object.text("id", "1ABC").number("atoms", 1_960);
    assert_eq!(object.finish(), r#"{"id":"1ABC","atoms":1960}"#);
}

#[test]
fn an_empty_object_renders_as_one() {
    assert_eq!(Json::new().finish(), "{}");
    assert_eq!(json_array(&[]), "[]");
}

#[test]
fn characters_a_string_cannot_carry_literally_are_escaped() {
    let mut object = Json::new();
    object.text("title", "a \"quoted\" \\ line\nbreak\ttab");
    assert_eq!(
        object.finish(),
        r#"{"title":"a \"quoted\" \\ line\nbreak\ttab"}"#
    );
}

#[test]
fn a_control_character_is_escaped_by_its_code_point() {
    let mut object = Json::new();
    object.text("odd", "\u{7}");
    assert_eq!(object.finish(), "{\"odd\":\"\\u0007\"}");
}

#[test]
fn a_list_of_objects_renders_as_an_array() {
    let mut first = Json::new();
    first.text("chain", "A");
    let mut second = Json::new();
    second.text("chain", "B");
    assert_eq!(
        json_array(&[first.finish(), second.finish()]),
        r#"[{"chain":"A"},{"chain":"B"}]"#
    );
}
