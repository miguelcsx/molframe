use super::*;

#[test]
fn the_bundled_reach_table_is_the_repository_copy() {
    let root = format!(
        "{}/../../data/macromolecular-reach.reach",
        env!("CARGO_MANIFEST_DIR")
    );
    let beside = format!(
        "{}/data/macromolecular-reach.reach",
        env!("CARGO_MANIFEST_DIR")
    );
    let (Ok(root), Ok(beside)) = (std::fs::read(&root), std::fs::read(&beside)) else {
        panic!("both copies are readable");
    };
    assert_eq!(root, beside);
}

#[cfg(feature = "embedded-tables")]
#[test]
fn the_bundled_profile_names_a_reach_table_that_decodes() {
    use crate::BondReach;
    let Some(tables) = bundled() else {
        panic!("the feature compiles the layer in");
    };
    let Some(name) = tables.profile_entry(<BondReach as crate::TableKind>::KIND) else {
        panic!("the profile names the reach");
    };
    let Some(entry) = tables.entry(&name, <BondReach as crate::TableKind>::EXTENSION) else {
        panic!("the named entry is held");
    };
    let Ok(provided) = entry.provided::<BondReach>() else {
        panic!("the reach decodes");
    };
    assert_eq!(&provided.table, crate::test_data::reach());
}

#[cfg(feature = "embedded-tables")]
#[test]
fn the_bundled_layer_is_built_once() {
    let (Some(first), Some(second)) = (bundled(), bundled()) else {
        panic!("the feature compiles the layer in");
    };
    assert!(Arc::ptr_eq(&first, &second));
}

#[cfg(not(feature = "embedded-tables"))]
#[test]
fn without_the_feature_nothing_is_compiled_in() {
    assert!(bundled().is_none());
}
