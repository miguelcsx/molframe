use super::*;

#[test]
fn a_spaced_pattern_projects_only_selected_positions() {
    let Ok(pattern) = SeedPattern::new(&[true, false, true, false]) else {
        panic!("nonempty pattern");
    };
    let seeds = pattern.seeds(b"ABCDE").collect::<Vec<_>>();
    assert_eq!(seeds, vec![(0, b"AC".to_vec()), (1, b"BD".to_vec())]);
}

#[test]
fn an_empty_weight_is_rejected() {
    assert_eq!(SeedPattern::new(&[false, false]), Err(SeedPatternError));
}
