use super::*;

fn members(selection: &AtomSelection) -> Vec<u32> {
    selection.iter().collect()
}

#[test]
fn a_contiguous_run_is_stored_as_a_run_rather_than_as_positions() {
    let selection = AtomSelection::from_sorted((10..20).collect());
    assert!(matches!(selection, AtomSelection::Range(_)));
    assert_eq!(selection.heap_bytes(), 0);
    assert_eq!(members(&selection), (10..20).collect::<Vec<_>>());
}

#[test]
fn a_few_runs_are_stored_as_runs() {
    let positions: Vec<u32> = (0..5).chain(100..105).collect();
    let selection = AtomSelection::from_sorted(positions.clone());
    assert!(matches!(selection, AtomSelection::Ranges(_)));
    assert_eq!(members(&selection), positions);
}

#[test]
fn a_scattered_handful_stays_a_list_rather_than_becoming_a_mask() {
    let selection = AtomSelection::from_sorted(vec![1, 5000, 999_999]);
    assert!(matches!(selection, AtomSelection::Sparse(_)));
    assert!(selection.heap_bytes() < 100);
}

#[test]
fn a_dense_scatter_becomes_a_mask() {
    let positions: Vec<u32> = (0..1000).filter(|value| value % 3 != 0).collect();
    let selection = AtomSelection::from_sorted(positions.clone());
    assert!(matches!(selection, AtomSelection::Dense(_)));
    assert_eq!(members(&selection), positions);
}

#[test]
fn intersection_of_two_runs_stays_a_run_without_enumerating_either() {
    let left = AtomSelection::range(0..100);
    let right = AtomSelection::range(50..200);
    let both = left.intersect(&right);
    assert_eq!(both, AtomSelection::Range(50..100));
}

#[test]
fn the_algebra_agrees_with_set_operations_on_every_shape() {
    let shapes = [
        AtomSelection::Empty,
        AtomSelection::All(64),
        AtomSelection::range(4..12),
        AtomSelection::from_sorted(vec![1, 2, 3, 40, 41, 60]),
        AtomSelection::from_sorted((0..64).filter(|value| value % 2 == 0).collect()),
        AtomSelection::from_sorted(vec![0, 63]),
    ];
    for left in &shapes {
        for right in &shapes {
            let (l, r) = (members(left), members(right));

            let expected_and: Vec<u32> = l
                .iter()
                .copied()
                .filter(|value| r.contains(value))
                .collect();
            assert_eq!(
                members(&left.intersect(right)),
                expected_and,
                "{left:?} & {right:?}"
            );

            let mut expected_or = l.clone();
            expected_or.extend(r.iter().copied());
            expected_or.sort_unstable();
            expected_or.dedup();
            assert_eq!(
                members(&left.union(right)),
                expected_or,
                "{left:?} | {right:?}"
            );

            let expected_not: Vec<u32> = l
                .iter()
                .copied()
                .filter(|value| !r.contains(value))
                .collect();
            assert_eq!(
                members(&left.difference(right)),
                expected_not,
                "{left:?} - {right:?}"
            );
        }
    }
}

#[test]
fn an_intersection_is_never_larger_than_its_smaller_operand() {
    let left = AtomSelection::range(0..1000);
    let right = AtomSelection::from_sorted(vec![10, 20, 30, 5000]);
    assert!(left.intersect(&right).len() <= right.len());
    assert!(left.intersect(&right).len() <= left.len());
}

#[test]
fn a_complement_covers_exactly_what_the_selection_does_not() {
    let selection = AtomSelection::from_sorted(vec![0, 2, 4]);
    let rest = selection.complement(6);
    assert_eq!(members(&rest), [1, 3, 5]);
    assert_eq!(selection.union(&rest).len(), 6);
    assert!(selection.intersect(&rest).is_empty());
}

#[test]
fn membership_agrees_with_enumeration_for_every_shape() {
    let selection = AtomSelection::from_sorted(vec![3, 4, 5, 90]);
    for position in 0..100u32 {
        assert_eq!(
            selection.contains(position),
            members(&selection).contains(&position),
            "position {position}"
        );
    }
}

#[test]
fn building_from_unsorted_positions_sorts_and_deduplicates() {
    let selection: AtomSelection = [5u32, 1, 5, 3, 1].into_iter().collect();
    assert_eq!(members(&selection), [1, 3, 5]);
}

#[test]
fn an_empty_run_collapses_rather_than_being_stored_as_a_run() {
    assert_eq!(AtomSelection::range(5..5), AtomSelection::Empty);
    let (high, low) = (9u32, 2u32);
    assert_eq!(AtomSelection::range(high..low), AtomSelection::Empty);
    assert!(AtomSelection::from_sorted(Vec::new()).is_empty());
}

#[test]
fn every_selection_shape_iterates_ascending_and_reports_its_exact_length() {
    let dense = {
        let mut mask = crate::column::BitVec::repeat(false, 200);
        for position in [3u32, 64, 65, 199] {
            mask.set(position, true);
        }
        AtomSelection::Dense(mask)
    };

    let cases = [
        (AtomSelection::Empty, Vec::new()),
        (AtomSelection::All(4), vec![0, 1, 2, 3]),
        (AtomSelection::Range(2..5), vec![2, 3, 4]),
        (
            AtomSelection::from_sorted(vec![1, 2, 3, 9, 10]),
            vec![1, 2, 3, 9, 10],
        ),
        (AtomSelection::from_sorted(vec![0, 5, 900]), vec![0, 5, 900]),
        (dense, vec![3, 64, 65, 199]),
    ];

    for (selection, expected) in cases {
        let iterator = selection.iter();
        let (lower, upper) = iterator.size_hint();
        assert_eq!(lower, expected.len(), "lower bound for {selection:?}");
        assert_eq!(upper, Some(expected.len()), "upper bound for {selection:?}");

        let collected: Vec<u32> = selection.iter().collect();
        assert_eq!(collected, expected, "positions of {selection:?}");
        assert!(
            collected.windows(2).all(|pair| pair[0] < pair[1]),
            "positions of {selection:?} must ascend"
        );

        // A fused iterator keeps returning nothing once it is finished.
        let mut drained = selection.iter();
        while drained.next().is_some() {}
        assert_eq!(drained.next(), None);
    }
}
