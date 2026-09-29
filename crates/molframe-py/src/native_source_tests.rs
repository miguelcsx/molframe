use super::{NativeBond, NativeStructureSource};
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

fn imported(py: Python<'_>) -> NativeStructureSource {
    let object = match Bound::new(py, PyStructure::new(structure())) {
        Ok(value) => value,
        Err(error) => panic!("structure should bind: {error}"),
    };
    match NativeStructureSource::from_python(object.as_any()) {
        Ok(value) => value,
        Err(error) => panic!("native source should import: {error}"),
    }
}

#[test]
fn a_sparse_selection_allocates_only_the_selected_rows() {
    Python::initialize();
    Python::attach(|py| {
        let source = imported(py);
        let rows = match source.select("index 1") {
            Ok(value) => value,
            Err(error) => panic!("query should evaluate: {error}"),
        };
        assert_eq!(rows, [1]);
        assert_eq!(rows.capacity(), 1);
        let none = match source.select("none") {
            Ok(value) => value,
            Err(error) => panic!("query should evaluate: {error}"),
        };
        assert!(none.is_empty());
        assert_eq!(none.capacity(), 0);
    });
}

#[test]
fn repeated_and_interleaved_queries_return_their_own_rows() {
    Python::initialize();
    Python::attach(|py| {
        let source = imported(py);
        for _ in 0..3 {
            for (query, expected) in [
                ("index 0", vec![0]),
                ("all", vec![0, 1]),
                ("index 1", vec![1]),
            ] {
                let rows = match source.select(query) {
                    Ok(value) => value,
                    Err(error) => panic!("{query} should evaluate: {error}"),
                };
                assert_eq!(rows, expected, "{query}");
            }
        }
    });
}

#[test]
fn a_malformed_query_is_a_value_error_and_leaves_the_source_usable() {
    Python::initialize();
    Python::attach(|py| {
        let source = imported(py);
        assert!(source.select("resname (").is_err());
        assert!(matches!(source.select("all"), Ok(rows) if rows == [0, 1]));
    });
}

#[test]
fn native_topology_preserves_file_bond_order() {
    const MMCIF: &str = "data_bond\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 N N GLY A 1 1 0 0\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n1 covale A 1 GLY CA A 1 GLY N SING\n";
    let input = molframe::InputBuffer::from_bytes(MMCIF.as_bytes().to_vec());
    let structure =
        match molframe::read_buffer(&input, Some("bonded.cif"), &molframe::ReadOptions::new()) {
            Ok((structure, _)) => structure,
            Err(findings) => panic!("bond fixture should parse: {findings:?}"),
        };
    Python::initialize();
    Python::attach(|py| {
        let object = match Bound::new(py, PyStructure::new(structure)) {
            Ok(value) => value,
            Err(error) => panic!("structure should bind: {error}"),
        };
        let source = match NativeStructureSource::from_python(object.as_any()) {
            Ok(value) => value,
            Err(error) => panic!("native source should import: {error}"),
        };
        let bonds = source.topology().bonds;
        assert_eq!(bonds.len(), 1);
        assert_eq!(bonds[0].order, NativeBond::ORDER_SINGLE);
        assert_eq!(bonds[0].aromatic, 0);
    });
}
