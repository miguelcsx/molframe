use super::*;

fn built(degrees: &[usize], items: &[(usize, u32)]) -> Csr<u32> {
    let mut builder = CsrBuilder::with_degrees(degrees);
    for (row, item) in items {
        builder.push(*row, *item);
    }
    builder.finish()
}

#[test]
fn an_empty_adjacency_has_no_rows_and_no_items() {
    let csr = Csr::<u32>::empty();

    assert!(csr.is_empty());
    assert_eq!(csr.rows(), 0);
    assert_eq!(csr.len(), 0);
    assert_eq!(csr.row(0), &[] as &[u32]);
}

#[test]
fn each_row_returns_exactly_the_items_scattered_into_it() {
    let csr = built(&[2, 0, 3], &[(0, 7), (2, 1), (0, 9), (2, 2), (2, 3)]);

    assert_eq!(csr.rows(), 3);
    assert_eq!(csr.len(), 5);
    assert_eq!(csr.row(0), &[7, 9]);
    assert_eq!(csr.row(1), &[] as &[u32]);
    assert_eq!(csr.row(2), &[1, 2, 3]);
}

#[test]
fn a_row_past_the_end_is_empty_rather_than_a_panic() {
    let csr = built(&[1], &[(0, 4)]);

    assert_eq!(csr.row(1), &[] as &[u32]);
    assert_eq!(csr.row(usize::MAX), &[] as &[u32]);
}

#[test]
fn items_beyond_a_declared_degree_never_reach_a_neighbouring_row() {
    let csr = built(&[1, 1], &[(0, 5), (0, 6), (0, 7), (1, 8)]);

    assert_eq!(csr.row(0), &[5]);
    assert_eq!(csr.row(1), &[8]);
}

#[test]
fn a_row_beyond_the_declared_count_is_dropped() {
    let csr = built(&[1], &[(0, 5), (9, 6)]);

    assert_eq!(csr.rows(), 1);
    assert_eq!(csr.row(0), &[5]);
}

#[test]
fn rows_can_be_sorted_in_place_through_row_mut() {
    let mut csr = built(&[3], &[(0, 9), (0, 2), (0, 5)]);
    csr.row_mut(0).sort_unstable();

    assert_eq!(csr.row(0), &[2, 5, 9]);
    assert_eq!(csr.row_mut(4), &mut [] as &mut [u32]);
}

#[test]
fn iterating_yields_every_row_in_order() {
    let csr = built(&[2, 0, 1], &[(0, 1), (0, 2), (2, 3)]);
    let rows: Vec<&[u32]> = csr.iter().collect();

    assert_eq!(rows, vec![&[1, 2][..], &[][..], &[3][..]]);
}

#[test]
fn a_degree_table_of_zeros_allocates_rows_but_no_items() {
    let csr = built(&[0, 0, 0], &[]);

    assert_eq!(csr.rows(), 3);
    assert_eq!(csr.len(), 0);
    assert!(csr.iter().all(<[u32]>::is_empty));
}

#[test]
fn an_under_filled_row_leaves_no_default_padding_behind() {
    let csr = built(&[3, 2], &[(0, 7), (1, 4)]);

    assert_eq!(csr.row(0), &[7]);
    assert_eq!(csr.row(1), &[4]);
    assert_eq!(csr.len(), 2);
}

#[test]
fn canonicalising_sorts_each_row_and_drops_repeated_items() {
    let mut builder = CsrBuilder::with_degrees(&[4, 2]);
    for (row, item) in [(0, 9), (0, 2), (0, 9), (0, 2), (1, 5), (1, 5)] {
        builder.push(row, item);
    }
    let csr = builder.finish_canonical();

    assert_eq!(csr.row(0), &[2, 9]);
    assert_eq!(csr.row(1), &[5]);
    assert_eq!(csr.len(), 3);
}

#[test]
fn canonicalising_keeps_distinct_items_that_repeat_across_rows() {
    let mut builder = CsrBuilder::with_degrees(&[2, 2]);
    for (row, item) in [(0, 1), (0, 1), (1, 1), (1, 3)] {
        builder.push(row, item);
    }
    let csr = builder.finish_canonical();

    assert_eq!(csr.row(0), &[1]);
    assert_eq!(csr.row(1), &[1, 3]);
}
