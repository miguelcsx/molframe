use super::*;
use std::path::PathBuf;

#[test]
fn compiled_query_delegates_selection_to_the_rust_facade() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match molframe::read(path) {
        Ok(structure) => PyStructure::new(structure),
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    };
    Python::initialize();
    Python::attach(|py| {
        let query = match PyQuery::new(py, "name CA") {
            Ok(query) => query,
            Err(error) => panic!("query compilation failed: {error}"),
        };
        let selected = match structure.select(py, PyQueryArg::Compiled(query), None, None) {
            Ok(selected) => selected,
            Err(error) => panic!("query evaluation failed: {error}"),
        };
        assert_eq!(selected.inner.iter().collect::<Vec<_>>(), vec![1]);
    });
}

#[test]
fn text_queries_compile_in_rust_and_reject_bad_syntax() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match molframe::read(path) {
        Ok(structure) => PyStructure::new(structure),
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    };
    Python::initialize();
    Python::attach(|py| {
        let selected = match structure.select(py, PyQueryArg::Text("name CA".into()), None, None) {
            Ok(selected) => selected,
            Err(error) => panic!("text query failed: {error}"),
        };
        assert_eq!(selected.inner.iter().collect::<Vec<_>>(), vec![1]);
        assert!(
            structure
                .select(py, PyQueryArg::Text("name".into()), None, None)
                .is_err()
        );
    });
}
