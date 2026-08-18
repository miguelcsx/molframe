use super::audit_batch;
use crate::{PolicyDimension, PolicySpace};
use pdbiox_core::contract::{AnalysisPolicy, MissingPolicy, PolicyField};
use std::collections::BTreeSet;

fn plan() -> crate::AuditPlan {
    match PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::missing_atoms([
            MissingPolicy::Report,
            MissingPolicy::Indeterminate,
        ]))
        .plan()
    {
        Ok(plan) => plan,
        Err(error) => panic!("plan failed: {error}"),
    }
}

#[test]
fn stability_is_pooled_across_subjects() {
    // Subject 0 is policy-independent (fully stable); subject 1 gains an item
    // only under the Report policy (stability 1/2).
    let subjects = [0u32, 1u32];
    let result = audit_batch(
        &plan(),
        &subjects,
        |subject: &u32, policy: &AnalysisPolicy| {
            let mut items = BTreeSet::from(["base"]);
            if *subject == 1 && policy.missing_atoms == MissingPolicy::Report {
                items.insert("extra");
            }
            Ok::<_, ()>(items)
        },
        Clone::clone,
    );
    let Ok(result) = result else {
        panic!("audit failed");
    };
    assert_eq!(result.subjects, 2);
    assert!(
        (result.mean_stability - 0.75).abs() < 1e-12,
        "{}",
        result.mean_stability
    );
    assert!((result.min_stability - 0.5).abs() < 1e-12);
    assert!((result.max_stability - 1.0).abs() < 1e-12);
    assert_eq!(result.dimensions.len(), 1);
    assert_eq!(result.dimensions[0].field, PolicyField::MissingAtoms);
    assert!(result.dimensions[0].mean_change > 0.0);
}

#[test]
fn an_empty_subject_set_is_fully_stable() {
    let empty: [u32; 0] = [];
    let result = audit_batch(
        &plan(),
        &empty,
        |_: &u32, _: &AnalysisPolicy| Ok::<_, ()>(BTreeSet::<&str>::new()),
        Clone::clone,
    );
    let Ok(result) = result else {
        panic!("audit failed");
    };
    assert_eq!(result.subjects, 0);
    assert!((result.mean_stability - 1.0).abs() < f64::EPSILON);
    assert!(result.dimensions.is_empty());
}
