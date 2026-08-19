use super::PyPolymerKind;

#[test]
fn polymer_kind_none_has_a_python_safe_constructor() {
    assert_eq!(PyPolymerKind::none(), PyPolymerKind::None);
}
