use super::{Comparison, MissingVerdict, VerdictProfile, VerdictRule, VerdictStatus};
use std::collections::BTreeMap;

fn profile(missing: MissingVerdict) -> VerdictProfile {
    VerdictProfile::new(
        "molframe-fx-strict-1.0",
        [
            VerdictRule {
                metric: "required_atom_rmsd".into(),
                comparison: Comparison::AtMost(1.0),
            },
            VerdictRule {
                metric: "coverage".into(),
                comparison: Comparison::AtLeast(0.95),
            },
        ],
        missing,
    )
}

#[test]
fn decisions_retain_the_exact_profile_and_each_rule() {
    let metrics = BTreeMap::from([
        (Box::<str>::from("coverage"), 1.0),
        (Box::<str>::from("required_atom_rmsd"), 0.8),
    ]);
    let verdict = profile(MissingVerdict::Indeterminate).decide(&metrics);
    assert_eq!(verdict.profile.as_ref(), "molframe-fx-strict-1.0");
    assert_eq!(verdict.status, VerdictStatus::Pass);
    assert_eq!(verdict.outcomes.len(), 2);
}

#[test]
fn absent_metrics_are_indeterminate_instead_of_invented() {
    let verdict = profile(MissingVerdict::Indeterminate).decide(&BTreeMap::new());
    assert_eq!(verdict.status, VerdictStatus::Indeterminate);
    assert!(
        verdict
            .outcomes
            .iter()
            .all(|outcome| outcome.value.is_none())
    );
}
