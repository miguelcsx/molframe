use super::*;

#[test]
fn a_dtype_mismatch_names_both_dtypes_and_the_conversion() {
    let message = dtype_mismatch("float32", "float64");
    assert!(message.contains("expected a float32 array"));
    assert!(message.contains("received dtype float64"));
    assert!(message.contains("astype(numpy.float32)"));
}

#[test]
fn a_rank_mismatch_names_both_ranks() {
    let message = rank_mismatch(2, 1);
    assert!(message.contains("2-dimensional"));
    assert!(message.contains("received 1 dimensions"));
}

#[test]
fn a_non_array_names_the_type_received() {
    assert!(not_an_array("list").contains("received list"));
}
