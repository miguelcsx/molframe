use super::*;

#[test]
fn percent_handles_empty_and_ratio() {
    assert!(percent(1, 0).is_nan());
    assert!((percent(1, 4) - 25.0).abs() < 1e-12);
}

#[test]
fn nucleic_counts_skip_unevaluated() {
    let table = NucleicTable {
        nu: vec![],
        pucker_phase: vec![],
        pucker_amplitude: vec![],
        pucker_class: vec![PUCKER_C3_ENDO, PUCKER_C2_ENDO, PUCKER_OTHER, NOT_EVALUATED],
        torsions: vec![],
        suite_flags: vec![],
        suite_class: vec![0, 1, NOT_EVALUATED, 1],
        tables: vec![],
    };
    let s = table.summary();
    assert_eq!(
        (s.pucker_evaluated, s.c3_endo, s.c2_endo, s.other_pucker),
        (3, 1, 1, 1)
    );
    assert_eq!((s.suite_evaluated, s.suite_flagged), (3, 2));
}
