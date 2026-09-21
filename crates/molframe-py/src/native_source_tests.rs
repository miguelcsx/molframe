use super::NativeStructureSource;
use crate::bindings::PyStructure;
use pyo3::prelude::*;

const PDB: &str = "\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
END
";

fn structure() -> molframe::Structure {
    let input = molframe::InputBuffer::from_bytes(PDB.as_bytes().to_vec());
    let result = molframe::read_buffer(&input, Some("source.pdb"), &molframe::ReadOptions::new());
    match result {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture should parse: {findings:?}"),
    }
}

#[test]
fn native_source_retains_coordinates_and_evaluates_queries() {
    Python::initialize();
    Python::attach(|py| {
        let structure = structure();
        let expected = structure.coordinates().as_ptr();
        let object = match Bound::new(py, PyStructure::new(structure)) {
            Ok(value) => value,
            Err(error) => panic!("structure should bind: {error}"),
        };
        let source = match NativeStructureSource::from_python(object.as_any()) {
            Ok(value) => value,
            Err(error) => panic!("native source should import: {error}"),
        };

        assert_eq!(source.coordinates().as_ptr(), expected);
        assert_eq!(source.coordinates().len(), 2);
        assert_eq!(source.topology().atoms.len(), 2);
        assert_eq!(source.topology().residue_atom_start, [0, 2]);
        let rows = match source.select("name CA") {
            Ok(value) => value,
            Err(error) => panic!("query should evaluate: {error}"),
        };
        assert_eq!(rows, [1]);
    });
}
