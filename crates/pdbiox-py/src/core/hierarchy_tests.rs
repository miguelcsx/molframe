use super::{PyChains, PyModels, PyResidues};
use pyo3::Python;
use std::path::PathBuf;

fn sample() -> pdbiox::Structure {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    match pdbiox::read(path) {
        Ok(structure) => structure,
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    }
}

#[test]
fn hierarchy_collections_delegate_to_rust_handles() {
    let structure = sample();
    let models = PyModels {
        inner: structure.clone(),
    };
    let chains = PyChains {
        inner: structure.clone(),
        model: None,
    };
    let residues = PyResidues {
        inner: structure,
        chain: None,
    };

    assert_eq!(models.__len__(), 1);
    assert_eq!(chains.__len__(), 1);
    assert_eq!(residues.__len__(), 1);
    Python::initialize();
    Python::attach(|py| {
        let chain = match chains.chain_named(py, "A") {
            Ok(chain) => chain,
            Err(error) => panic!("chain lookup failed: {error}"),
        };
        assert_eq!(chain.label().as_deref(), Some("A"));
    });
}
