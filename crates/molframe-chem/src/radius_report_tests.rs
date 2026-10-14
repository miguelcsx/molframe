use super::*;
use crate::test_data::{bondi, metal_tables};

#[test]
fn the_policy_names_every_table_behind_it() {
    let policy = RadiusPolicy {
        source: bondi(),
        metals: MetalFallback::Ionic,
        metal_tables: Some(metal_tables()),
        hetero: HeteroPolicy::POLYMER,
    };
    let applied = policy.applied();
    assert_eq!(applied.radii.name, "bondi");
    assert_eq!(applied.metal_tables.len(), 2);
    let text = applied.to_string();
    assert!(text.contains("bondi 1964"), "{text}");
    assert!(text.contains("shannon-ionic 1976"), "{text}");
    assert!(text.contains("hetero none"), "{text}");
}

#[test]
fn no_fallback_names_no_metal_table() {
    let policy = RadiusPolicy {
        source: bondi(),
        metals: MetalFallback::None,
        metal_tables: None,
        hetero: HeteroPolicy::POLYMER,
    };
    assert!(policy.applied().metal_tables.is_empty());
}
