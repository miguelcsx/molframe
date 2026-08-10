use super::*;

#[test]
fn stars_questions_classes_negation_ranges_and_escaping_match() {
    assert!(Glob::new("C*").matches("CA"));
    assert!(Glob::new("C?").matches("CB"));
    assert!(Glob::new("[CN]A").matches("CA"));
    assert!(!Glob::new("[!CN]A").matches("CA"));
    assert!(Glob::new("[A-C]A").matches("BA"));
    assert!(Glob::new("C\\*").matches("C*"));
}

#[test]
fn a_star_can_match_empty_and_several_characters() {
    let glob = Glob::new("A*B");
    assert!(glob.matches("AB"));
    assert!(glob.matches("ALONGB"));
    assert!(!glob.matches("ALONG"));
}
