use super::*;

#[test]
fn every_named_role_is_unique_and_decodes_back() {
    let mut codes = std::collections::BTreeSet::new();
    for (name, role) in PolymerAtomRole::NAMED {
        assert!(codes.insert(role.code()), "{name} repeats a code");
        assert_eq!(
            PolymerAtomRole::from_code(role.code()),
            Some(role),
            "{name}"
        );
    }
}
