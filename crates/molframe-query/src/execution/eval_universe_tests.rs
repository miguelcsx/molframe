use super::tests::structure;
use super::*;

#[test]
fn a_universe_reaching_past_the_structure_is_rejected_before_anything_is_selected() {
    let structure = structure();
    let query = Query::compile("all").expect("query compiles");
    let policy = AnalysisPolicy::default();
    let past_end = AtomSelection::Sparse(vec![0, structure.atom_count()]);
    let result = query.evaluate_in(&structure, &past_end, &policy, &Groups::new(), None);
    match result {
        Err(findings) => assert!(
            findings.iter().any(|finding| finding.code() == Code::E6009),
            "an out-of-range universe reports the out-of-bounds atom: {findings:?}"
        ),
        Ok(_) => panic!("a universe naming an atom past the structure must be rejected"),
    }
}

#[test]
fn a_universe_covering_the_whole_structure_is_not_read_as_out_of_range() {
    let structure = structure();
    let query = Query::compile("all").expect("query compiles");
    let policy = AnalysisPolicy::default();
    let all = AtomSelection::All(structure.atom_count());
    let result = query
        .evaluate_in(&structure, &all, &policy, &Groups::new(), None)
        .unwrap_or_else(|findings| panic!("a covering universe must evaluate: {findings:?}"));
    assert_eq!(
        result.selection.iter().collect::<Vec<_>>(),
        (0..structure.atom_count()).collect::<Vec<_>>(),
        "every atom is selected"
    );
}

#[test]
fn an_empty_universe_is_a_valid_universe() {
    let structure = structure();
    let query = Query::compile("all").expect("query compiles");
    let policy = AnalysisPolicy::default();
    let result = query
        .evaluate_in(
            &structure,
            &AtomSelection::Empty,
            &policy,
            &Groups::new(),
            None,
        )
        .unwrap_or_else(|findings| panic!("an empty universe must evaluate: {findings:?}"));
    assert!(
        result.selection.is_empty(),
        "nothing is inside an empty universe"
    );
}
