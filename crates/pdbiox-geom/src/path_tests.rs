use super::path_torsions;

#[test]
fn missing_atoms_invalidate_only_their_dependent_path_torsions() {
    let atoms = [
        Some([0.0, 0.0, 0.0]),
        Some([1.0, 0.0, 0.0]),
        Some([1.0, 1.0, 0.0]),
        Some([1.0, 1.0, 1.0]),
        None,
        Some([2.0, 2.0, 2.0]),
        Some([3.0, 2.0, 2.0]),
    ];
    let torsions = path_torsions(&atoms, 5);
    assert_eq!(torsions.len(), 4);
    assert!(torsions[0].is_some());
    assert!(torsions[1..].iter().all(Option::is_none));
}
