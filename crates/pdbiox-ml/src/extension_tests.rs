use super::*;

#[test]
fn registry_defines_every_documented_domain_extension() {
    let extensions = [
        PdbioxExtension::AtomIndex,
        PdbioxExtension::ResidueIndex,
        PdbioxExtension::ChainIndex,
        PdbioxExtension::EntityIndex,
        PdbioxExtension::Coordinates3f,
        PdbioxExtension::Element,
        PdbioxExtension::SymbolId,
        PdbioxExtension::Altloc,
        PdbioxExtension::Selection,
        PdbioxExtension::Validity,
    ];
    assert_eq!(extensions.len(), 10);
    assert!(extensions.iter().all(|extension| {
        extension.name().starts_with("pdbiox.") && extension.storage_type() != DataType::Null
    }));
}
