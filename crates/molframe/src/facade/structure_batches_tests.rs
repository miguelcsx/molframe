use super::*;
use molframe_core::{ScratchPolicy, TempStoragePolicy, collect_structure};
use std::io::Write;

const PDB: &str = concat!(
    "ATOM      1  N   GLY A   1      11.104  13.207  10.111  1.00 20.00           N  \n",
    "ATOM      2  CA  GLY A   1      12.104  13.207  10.111  1.00 21.00           C  \n",
    "END\n",
);

fn context() -> ExecutionContext {
    ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context")
}

#[test]
#[cfg(feature = "pdb")]
fn dispatch_and_collector_never_reserve_from_the_input_length() {
    let mut file = tempfile::Builder::new()
        .suffix(".pdb")
        .tempfile()
        .expect("temporary input");
    file.write_all(PDB.as_bytes()).expect("fixture write");
    let context = context();
    let mut source = open_structure_batches(file.path(), &ReadOptions::new(), &context)
        .expect("bounded PDB source");
    let (structure, diagnostics) =
        collect_structure(&mut source, &context).expect("budgeted collector");

    assert_eq!(structure.atom_count(), 2);
    assert!(diagnostics.is_empty());
    assert!(context.peak_reserved_bytes() < 16 * 1024 * 1024);
    assert_eq!(context.live_batches(), 0);
}

#[test]
#[cfg(all(feature = "bcif", feature = "pdb"))]
fn binary_cif_dispatch_uses_the_bounded_reader() {
    let input = molframe_core::InputBuffer::from_bytes(PDB.as_bytes().to_vec());
    let (structure, diagnostics) =
        molframe_pdb::read(&input, &ReadOptions::new()).expect("PDB structure fixture");
    assert!(diagnostics.is_empty());
    let write_options = molframe_cif::CifWriteOptions::new().with_block_id("batch");
    let bytes = molframe_bcif::write_structure_with_options(&structure, &write_options)
        .expect("BinaryCIF fixture");
    let mut file = tempfile::Builder::new()
        .suffix(".bcif")
        .tempfile()
        .expect("temporary BinaryCIF input");
    file.write_all(&bytes).expect("BinaryCIF fixture write");
    let context = context();
    let mut source = open_structure_batches(file.path(), &ReadOptions::new(), &context)
        .expect("bounded BinaryCIF source");
    let (decoded, findings) =
        collect_structure(&mut source, &context).expect("budgeted BinaryCIF collector");

    assert_eq!(decoded.atom_count(), structure.atom_count());
    assert!(findings.is_empty());
    assert!(context.peak_reserved_bytes() < 16 * 1024 * 1024);
    assert_eq!(context.live_batches(), 0);
}

#[test]
#[cfg(all(feature = "gzip", feature = "pdb"))]
fn compressed_batches_require_and_use_explicit_bounded_spill() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let mut file = tempfile::Builder::new()
        .suffix(".pdb.gz")
        .tempfile_in(directory.path())
        .expect("temporary gzip input");
    let mut encoder =
        flate2::write::GzEncoder::new(file.as_file_mut(), flate2::Compression::fast());
    encoder.write_all(PDB.as_bytes()).expect("gzip fixture");
    encoder.finish().expect("finish gzip fixture");

    let disabled = context();
    assert!(open_structure_batches(file.path(), &ReadOptions::new(), &disabled).is_err());

    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .temp_storage_policy(TempStoragePolicy::directory(
            directory.path().join("spill"),
            1024 * 1024,
        ))
        .build()
        .expect("spill context");
    let mut source = open_structure_batches(file.path(), &ReadOptions::new(), &context)
        .expect("bounded compressed source");
    assert!(context.spill_bytes() > 0);
    let (structure, diagnostics) =
        collect_structure(&mut source, &context).expect("compressed collection");
    assert_eq!(structure.atom_count(), 2);
    assert!(diagnostics.is_empty());
    drop(source);
    assert_eq!(context.spill_bytes(), 0);
}
