use super::*;
use crate::structure::fixture;

#[test]
fn a_view_of_everything_covers_every_atom() {
    let structure = fixture::sample();
    let view = structure.view();
    assert_eq!(view.len(), structure.atom_count());
    assert_eq!(view.len(), 24);
    assert!(!view.is_empty());
}

#[test]
fn narrowing_a_view_composes_without_touching_the_atoms_it_drops() {
    let structure = fixture::sample();
    let first_chain = structure.view_of(AtomSelection::range(0..12));
    let backbone = first_chain.narrow(&AtomSelection::from_sorted(
        (0..24).filter(|atom| atom % 4 == 1).collect(),
    ));
    assert_eq!(backbone.len(), 3, "one alpha carbon per residue of chain A");
}

#[test]
fn views_compose_as_sets() {
    let structure = fixture::sample();
    let left = structure.view_of(AtomSelection::range(0..12));
    let right = structure.view_of(AtomSelection::range(8..24));

    assert_eq!(left.union(&right).len(), 24);
    assert_eq!(left.intersect(&right).len(), 4);
    assert_eq!(left.difference(&right).len(), 8);
}

#[test]
fn a_view_yields_the_positions_of_the_atoms_it_covers() {
    let structure = fixture::sample();
    let view = structure.view_of(AtomSelection::from_sorted(vec![0, 5, 23]));
    let positions: Vec<_> = view.positions(ModelIndex::new(0)).collect();
    assert_eq!(positions.len(), 3);
    assert!(!view.bounds(ModelIndex::new(0)).is_empty());
}

#[test]
fn a_view_of_a_model_that_does_not_exist_yields_no_positions() {
    let structure = fixture::sample();
    let view = structure.view();
    assert_eq!(view.positions(ModelIndex::new(9)).count(), 0);
    assert!(view.bounds(ModelIndex::new(9)).is_empty());
}

#[test]
fn a_view_knows_whether_the_structure_has_moved_on_since_it_was_made() {
    let structure = fixture::sample();
    let view = structure.view();
    assert!(!view.is_stale_for(&structure));

    let mut moved = structure.data().clone();
    moved.generation = structure.generation().next();
    assert!(view.is_stale_for(&Structure::new(moved)));
}
