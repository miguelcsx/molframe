use super::{compact_ranges, resid_query};

#[test]
fn consecutive_numbers_collapse_into_ranges_whatever_the_order() {
    assert_eq!(compact_ranges(&[20, 21, 22, 25]), "20-22,25");
    assert_eq!(compact_ranges(&[5, 3, 4, 3]), "3-5");
    assert_eq!(compact_ranges(&[7]), "7");
    assert_eq!(compact_ranges(&[]), "");
    assert_eq!(compact_ranges(&[-2, -1, 0, 1]), "-2-1");
    assert_eq!(compact_ranges(&[i32::MAX - 1, i32::MAX]), "2147483646-2147483647");
}

#[test]
fn a_resid_query_joins_consecutive_numbers_into_ranges_per_chain() {
    let pairs: Vec<(String, i32)> = [("A", 3), ("A", 1), ("A", 2), ("A", 7), ("B", 4), ("A", 2)]
        .iter()
        .map(|(c, n)| ((*c).to_owned(), *n))
        .collect();
    let (query, count) = resid_query(&pairs).expect("residues were given");
    assert_eq!(count, 5);
    assert_eq!(
        query,
        "(chain A and (resid 1:3 or resid 7)) or (chain B and (resid 4))"
    );
}

#[test]
fn a_resid_query_for_no_residues_is_nothing_and_an_unnamed_chain_has_no_scope() {
    assert!(resid_query(&[]).is_none());
    let (query, count) = resid_query(&[(String::new(), 9)]).expect("one residue");
    assert_eq!((query.as_str(), count), ("(resid 9)", 1));
}
