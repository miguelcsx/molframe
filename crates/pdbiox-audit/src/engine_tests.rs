use super::audit;
use crate::{PolicyDimension, PolicySpace};
use pdbiox_core::contract::{AnalysisPolicy, MissingPolicy, Namespace, PolicyField};
use std::collections::BTreeSet;

#[test]
fn audit_identifies_global_and_dimension_specific_fragility() {
    let plan = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::identifiers([
            Namespace::Auth,
            Namespace::Label,
        ]))
        .vary(PolicyDimension::missing_atoms([
            MissingPolicy::Report,
            MissingPolicy::Indeterminate,
        ]))
        .plan()
        .unwrap_or_else(|error| panic!("plan failed: {error}"));
    let result = audit(
        &plan,
        |policy| {
            let mut items = BTreeSet::from(["stable"]);
            if policy.identifiers == Namespace::Auth {
                items.insert("namespace");
            }
            if policy.missing_atoms == MissingPolicy::Report {
                items.insert("coverage");
            }
            Ok::<_, ()>(items)
        },
        Clone::clone,
    )
    .unwrap_or_else(|()| panic!("analysis failed"));

    assert_eq!(result.runs.len(), 4);
    assert!((result.stability - 1.0 / 3.0).abs() < 1e-12);
    assert_eq!(result.sensitive_items.len(), 2);
    assert_eq!(result.dimensions[0].field, PolicyField::Identifiers);
    assert_eq!(result.dimensions[0].sensitive_items, ["namespace"]);
    assert_eq!(result.dimensions[1].field, PolicyField::MissingAtoms);
    assert_eq!(result.dimensions[1].sensitive_items, ["coverage"]);
    assert!((result.dimensions[0].mean_change - 5.0 / 12.0).abs() < 1e-12);
}

#[test]
fn empty_results_are_fully_stable() {
    let plan = PolicySpace::new(AnalysisPolicy::default())
        .plan()
        .unwrap_or_else(|error| panic!("plan failed: {error}"));
    let report = audit(&plan, |_| Ok::<_, ()>(BTreeSet::<u32>::new()), Clone::clone)
        .unwrap_or_else(|()| panic!("analysis failed"));
    assert!((report.stability - 1.0).abs() < f64::EPSILON);
    assert!(report.sensitive_items.is_empty());
}
