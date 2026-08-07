use super::normalise_index;

#[test]
fn python_indices_are_only_converted_at_the_boundary() {
    assert_eq!(normalise_index(0, 3), Some(0));
    assert_eq!(normalise_index(-1, 3), Some(2));
    assert_eq!(normalise_index(3, 3), None);
    assert_eq!(normalise_index(-4, 3), None);
}
