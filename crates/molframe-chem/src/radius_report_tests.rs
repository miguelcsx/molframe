use super::*;
use crate::{CustomRadii, RadiusSource};

#[test]
fn the_default_policy_names_every_table_behind_it() {
    let applied = RadiusPolicy::default().applied();
    assert_eq!(applied.radii.name, "bondi");
    assert_eq!(applied.metal_tables.len(), 2);
    let text = applied.to_string();
    assert!(text.contains("bondi 1964"), "{text}");
    assert!(text.contains("alvarez 2013"), "{text}");
    assert!(text.contains("hetero none"), "{text}");
}

#[test]
fn a_supplied_table_replaces_the_default_identity() {
    let Ok(table) = CustomRadii::parse("name: mine\nversion: 9\nsource: me\nlicence: MIT\nC 2.0\n")
    else {
        panic!("table parses");
    };
    let policy = RadiusPolicy {
        source: RadiusSource::from(table),
        metals: MetalFallback::None,
        ..RadiusPolicy::default()
    };
    let applied = policy.applied();
    assert_eq!(applied.radii.name, "mine");
    assert!(applied.metal_tables.is_empty());
}
