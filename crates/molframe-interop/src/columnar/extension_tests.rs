use super::*;

#[test]
fn registry_defines_every_documented_domain_extension() {
    let extensions = [
        MolframeExtension::AtomIndex,
        MolframeExtension::ResidueIndex,
        MolframeExtension::ChainIndex,
        MolframeExtension::EntityIndex,
        MolframeExtension::Coordinates3f,
        MolframeExtension::Element,
        MolframeExtension::SymbolId,
        MolframeExtension::Altloc,
        MolframeExtension::Selection,
        MolframeExtension::Validity,
    ];
    assert_eq!(extensions.len(), 10);
    assert!(extensions.iter().all(|extension| {
        extension.name().starts_with("molframe.") && extension.storage_type() != DataType::Null
    }));
}
