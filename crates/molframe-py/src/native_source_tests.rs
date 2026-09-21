use super::NativeStructureSource;
use crate::bindings::PyStructure;
use pyo3::prelude::*;

const MMCIF: &str = "\
data_source
_entry.id source
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N ALA A 1 1 11.104 6.134 -6.504 1.00 0.00
ATOM 2 C CA ALA A 1 1 12.560 6.195 -6.504 1.00 0.00
";

fn structure() -> molframe::Structure {
    let input = molframe::InputBuffer::from_bytes(MMCIF.as_bytes().to_vec());
    let result = molframe::read_buffer(&input, Some("source.cif"), &molframe::ReadOptions::new());
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
        let rows = match source.select("all") {
            Ok(value) => value,
            Err(error) => panic!("query should evaluate: {error}"),
        };
        assert_eq!(rows, [0, 1]);
        let first = match source.encode_bcif() {
            Ok(value) => value,
            Err(error) => panic!("BCIF should encode: {error}"),
        };
        let second = match source.encode_bcif() {
            Ok(value) => value,
            Err(error) => panic!("cached BCIF should encode: {error}"),
        };
        assert!(!first.is_empty());
        assert_eq!(first, second);
    });
}
