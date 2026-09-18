use super::{PlanError, PolicyDimension, PolicySpace};
use molframe_core::contract::{AnalysisPolicy, MissingPolicy, ModelChoice, Namespace};

#[test]
fn cost_is_available_before_expansion() {
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::model([
            ModelChoice::First,
            ModelChoice::All,
        ]))
        .vary(PolicyDimension::identifiers([
            Namespace::Auth,
            Namespace::Label,
        ]))
        .vary(PolicyDimension::missing_atoms([
            MissingPolicy::Report,
            MissingPolicy::Indeterminate,
        ]));
    assert_eq!(space.cost(), Ok(8));
    let Ok(plan) = space.plan() else {
        panic!("bounded space should plan");
    };
    assert_eq!(plan.cost(), 8);
    assert_eq!(plan.policies()[0].model, ModelChoice::First);
    assert_eq!(plan.policies()[7].model, ModelChoice::All);
}

#[test]
fn empty_duplicate_and_oversized_spaces_are_refused() {
    let empty = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::model([]))
        .cost();
    assert_eq!(
        empty,
        Err(PlanError::EmptyDimension(
            molframe_core::contract::PolicyField::Model
        ))
    );

    let duplicate = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::model([ModelChoice::First]))
        .vary(PolicyDimension::model([ModelChoice::All]))
        .cost();
    assert_eq!(
        duplicate,
        Err(PlanError::DuplicateDimension(
            molframe_core::contract::PolicyField::Model
        ))
    );

    let oversized = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::model([
            ModelChoice::First,
            ModelChoice::All,
        ]))
        .with_max_runs(1)
        .plan();
    assert!(matches!(
        oversized,
        Err(PlanError::LimitExceeded { cost: 2, limit: 1 })
    ));
}
