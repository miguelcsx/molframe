use super::{ReadOptions, read_bytes};

#[cfg(all(feature = "pdb", feature = "gzip"))]
#[test]
fn owned_compressed_bytes_are_decoded_before_format_detection() {
    let bytes = vec![
        31, 139, 8, 0, 0, 0, 0, 0, 2, 255, 115, 12, 241, 247, 85, 0, 3, 67, 5, 5, 103, 71, 5, 5,
        71, 31, 71, 5, 71, 8, 23, 12, 12, 244, 12, 12, 12, 48, 105, 67, 32, 165, 96, 9, 226, 40,
        32, 128, 179, 130, 2, 151, 171, 159, 11, 23, 0, 141, 81, 192, 56, 85, 0, 0, 0,
    ];
    let options = ReadOptions::new().only_atomic_coords(true);
    let (structure, _) = read_bytes(bytes.clone(), Some("model.pdb.gz"), &options).unwrap();
    assert_eq!(structure.atom_count(), 1);
    let context = molframe_core::ExecutionContext::default();
    let (budgeted, _) =
        super::read_bytes_in(bytes, Some("model.pdb.gz"), &options, &context).unwrap();
    assert_eq!(budgeted.atom_count(), structure.atom_count());
}
