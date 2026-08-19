use super::load_structure;
use pyo3::Python;
use std::path::PathBuf;

#[test]
fn compatibility_loading_delegates_to_the_native_facade() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let expected = match pdbiox::read(&path) {
        Ok(structure) => structure,
        Err(findings) => panic!("native fixture failed to read: {findings:?}"),
    };
    Python::initialize();
    Python::attach(|py| {
        let loaded = load_structure(py, path).expect("compatibility load should use native reader");
        assert_eq!(loaded.structure().atom_count(), expected.atom_count());
        assert_eq!(loaded.structure().positions(), expected.positions());
    });
}
