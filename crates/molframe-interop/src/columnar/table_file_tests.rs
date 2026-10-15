use crate::columnar::ipc::{write_atom_ipc, write_atom_ipc_with_metadata};
use crate::columnar::parquet::write_atom_parquet;
use crate::test_support::two_atoms as structure;
use arrow::ipc::reader::FileReader;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::collections::BTreeMap;
use std::fs::File;

#[test]
fn arrow_ipc_round_trip_preserves_atom_rows() {
    let destination = match tempfile::NamedTempFile::new() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = write_atom_ipc(destination.path(), &structure()) {
        panic!("IPC write failed: {error}")
    }
    let file = match File::open(destination.path()) {
        Ok(file) => file,
        Err(error) => panic!("IPC open failed: {error}"),
    };
    let reader = match FileReader::try_new(file, None) {
        Ok(reader) => reader,
        Err(error) => panic!("IPC read failed: {error}"),
    };
    let rows: usize = reader
        .map(|batch| batch.map_or(0, |batch| batch.num_rows()))
        .sum();
    assert_eq!(rows, 2);
}

#[test]
fn parquet_round_trip_preserves_atom_rows() {
    let destination = match tempfile::NamedTempFile::new() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = write_atom_parquet(destination.path(), &structure()) {
        panic!("Parquet write failed: {error}")
    }
    let file = match File::open(destination.path()) {
        Ok(file) => file,
        Err(error) => panic!("Parquet open failed: {error}"),
    };
    let builder = match ParquetRecordBatchReaderBuilder::try_new(file) {
        Ok(builder) => builder,
        Err(error) => panic!("Parquet reader failed: {error}"),
    };
    let reader = match builder.build() {
        Ok(reader) => reader,
        Err(error) => panic!("Parquet batch reader failed: {error}"),
    };
    let rows: usize = reader
        .map(|batch| batch.map_or(0, |batch| batch.num_rows()))
        .sum();
    assert_eq!(rows, 2);
}

#[test]
fn arrow_file_carries_caller_provenance_metadata() {
    let destination = match tempfile::NamedTempFile::new() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = write_atom_ipc_with_metadata(
        destination.path(),
        &structure(),
        BTreeMap::from([("molframe.policy".to_owned(), "explicit".to_owned())]),
    ) {
        panic!("IPC write failed: {error}")
    }
    let file = match File::open(destination.path()) {
        Ok(file) => file,
        Err(error) => panic!("IPC open failed: {error}"),
    };
    let reader = match FileReader::try_new(file, None) {
        Ok(reader) => reader,
        Err(error) => panic!("IPC read failed: {error}"),
    };
    assert_eq!(
        reader
            .schema()
            .metadata()
            .get("molframe.policy")
            .map(String::as_str),
        Some("explicit")
    );
}
