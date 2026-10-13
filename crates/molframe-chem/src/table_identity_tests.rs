use super::*;

#[test]
fn display_names_all_four_facts() {
    let identity = TableIdentity::builtin("bondi", "1964", "Bondi", "MIT");
    assert_eq!(
        identity.to_string(),
        "bondi 1964 (source: Bondi; licence: MIT)"
    );
}

#[test]
fn builtin_and_supplied_identities_with_the_same_text_are_equal() {
    let built = TableIdentity::builtin("a", "1", "s", "l");
    assert_eq!(built, TableIdentity::supplied("a", "1", "s", "l"));
}
