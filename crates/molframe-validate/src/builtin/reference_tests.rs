use super::*;
use crate::builtin_validation::{BuiltinError, bundled_reference};

#[test]
fn overlay_replaces_only_what_is_given() -> Result<(), BuiltinError> {
    let base = bundled_reference()?.clone();
    let info = TableInfo::new("mine", "2", "test", "MIT");
    let overrides = ReferenceOverrides::default()
        .with_geometry_text("bond\tALA\tN-CA\t1.5\t0.02\tX\n", info.clone());
    let set = base.clone().overlaid(overrides);
    assert_eq!(set.geometry.info(), &info);
    assert_eq!(set.geometry.len(), 1);
    assert_eq!(set.rama.info(), base.rama.info());
    assert_eq!(set.infos()[2], info);
    Ok(())
}

#[test]
fn table_info_displays_every_field() {
    let text = TableInfo::new("n", "v", "s", "l").to_string();
    assert_eq!(text, "n v (s; l)");
}
