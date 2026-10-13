use super::*;

#[test]
fn bundled_set_declares_source_version_and_licence() -> Result<(), BuiltinError> {
    let infos = bundled_reference()?.infos();
    assert_eq!(infos.len(), 4);
    for info in &infos {
        assert!(!info.name.is_empty());
        assert!(!info.version.is_empty());
        assert!(!info.source.is_empty());
        assert!(!info.licence.is_empty());
    }
    let rama = &infos[0];
    assert_eq!(rama.licence, "CC-BY-4.0");
    assert!(rama.source.contains("rlabduke/reference_data"));
    assert!(infos[2].licence.contains("LGPL-3.0"));
    Ok(())
}

#[test]
fn bundled_set_is_decoded_once() -> Result<(), BuiltinError> {
    assert!(std::ptr::eq(bundled_reference()?, bundled_reference()?));
    Ok(())
}
