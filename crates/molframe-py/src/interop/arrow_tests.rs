use super::*;
use arrow::ffi_stream::{ArrowArrayStreamReader, FFI_ArrowArrayStream};
use molframe::Structure;
use pyo3::types::PyCapsuleMethods;
use std::path::PathBuf;

fn structure() -> Structure {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    match molframe::read(path) {
        Ok(structure) => structure,
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    }
}

#[test]
fn atoms_publish_a_named_one_shot_arrow_stream_capsule() {
    Python::initialize();
    Python::attach(|py| {
        let atoms = PyAtoms::new(structure());
        let capsule = match atoms.__arrow_c_stream__(py, None) {
            Ok(capsule) => capsule,
            Err(error) => panic!("capsule export failed: {error}"),
        };
        let pointer = match capsule.pointer_checked(Some(c"arrow_array_stream")) {
            Ok(pointer) => pointer.cast::<FFI_ArrowArrayStream>(),
            Err(error) => panic!("named stream pointer absent: {error}"),
        };
        // SAFETY: `pointer_checked` verified the protocol name and returns the
        // writable stream value owned by this live capsule. `from_raw` moves it
        // out and nulls the capsule's release callback, preventing double drop.
        let mut reader = match unsafe { ArrowArrayStreamReader::from_raw(pointer.as_ptr()) } {
            Ok(reader) => reader,
            Err(error) => panic!("C stream import failed: {error}"),
        };
        let batch = match reader.next() {
            Some(Ok(batch)) => batch,
            Some(Err(error)) => panic!("C stream batch failed: {error}"),
            None => panic!("C stream returned no atom batch"),
        };
        assert_eq!(batch.num_rows(), atoms.structure().atom_count() as usize);
    });
}

#[test]
fn every_python_table_delegates_to_its_rust_arrow_adapter() {
    Python::initialize();
    Python::attach(|py| {
        let structure = structure();
        let residues = PyResidues::new_full(structure.clone());
        let chains = PyChains::new_full(structure.clone());
        let bonds = PyBonds::new(structure);
        assert!(residues.__arrow_c_stream__(py, None).is_ok());
        assert!(chains.__arrow_c_stream__(py, None).is_ok());
        assert!(bonds.__arrow_c_stream__(py, None).is_ok());
    });
}
