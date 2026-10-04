use super::{PlanError, PolicyDimension, PolicySpace};
use molframe_core::contract::{AnalysisPolicy, MissingPolicy, ModelChoice, Namespace, PolicyField};

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

#[test]
fn alternatives_written_as_words_mean_what_a_written_policy_means() {
    let Ok(dimension) = PolicyDimension::named(
        PolicyField::Assembly,
        &["asymmetric-unit", "biological:1", "crystal:8"],
    ) else {
        panic!("a valid vocabulary");
    };
    assert_eq!(dimension.len(), 3);
    let Ok(tolerances) = PolicyDimension::named(PolicyField::FloatTolerance, &["1e-6,1e-9"]) else {
        panic!("a valid tolerance");
    };
    assert_eq!(tolerances.len(), 1);
    assert!(PolicyDimension::named(PolicyField::Altloc, &["first", "nonsense"]).is_err());
    assert!(PolicyDimension::named(PolicyField::FloatTolerance, &["1e-6"]).is_err());
    assert!(PolicyDimension::named(PolicyField::FloatTolerance, &["-1,0"]).is_err());
}

#[test]
fn a_constrained_plan_drops_what_cannot_run_and_says_it_is_not_the_whole_product() {
    use molframe_core::contract::{AssemblyChoice, SymmetryPolicy};
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::assembly([
            AssemblyChoice::AsymmetricUnit,
            AssemblyChoice::Biological("1".into()),
        ]))
        .vary(PolicyDimension::symmetry([
            SymmetryPolicy::None,
            SymmetryPolicy::Crystallographic,
        ]));
    let Ok(plan) = space.clone().plan_constrained() else {
        panic!("two combinations survive");
    };
    assert_eq!(plan.cost(), 2);
    assert_eq!(plan.skipped(), 2);
    assert!(!plan.balanced());
    // A strict plan over the same space refuses instead.
    assert!(matches!(space.plan(), Err(PlanError::Conflict { .. })));

    let all_removed =
        PolicySpace::new(AnalysisPolicy::default()).vary(PolicyDimension::symmetry([
            SymmetryPolicy::Crystallographic,
        ]));
    assert!(matches!(
        all_removed.plan_constrained(),
        Err(PlanError::NoUniverse)
    ));
}

#[test]
fn a_forbidden_pair_is_removed_from_a_constrained_plan_and_refuses_a_strict_one() {
    use crate::PolicyValue;
    use molframe_core::contract::{AltlocPolicy, HydrogenPolicy};
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::hydrogens([
            HydrogenPolicy::ExplicitOnly,
            HydrogenPolicy::Exclude,
        ]))
        .vary(PolicyDimension::altloc([
            AltlocPolicy::ConformerConsistent,
            AltlocPolicy::First,
        ]))
        .forbid(
            PolicyValue::Hydrogens(HydrogenPolicy::Exclude),
            PolicyValue::Altloc(AltlocPolicy::First),
        );
    let Ok(plan) = space.clone().plan_constrained() else {
        panic!("three combinations survive");
    };
    assert_eq!(plan.cost(), 3);
    assert!(!plan.balanced());
    assert!(matches!(space.plan(), Err(PlanError::Conflict { run: 3 })));
}

#[test]
fn a_dimension_carries_its_class_and_the_case_for_varying_it() {
    use molframe_core::contract::HydrogenPolicy;
    let dimension =
        PolicyDimension::hydrogens([HydrogenPolicy::ExplicitOnly, HydrogenPolicy::Exclude])
            .justified(
                "structures differ in whether hydrogens were modelled",
                "PDBbind v2020 protein files carry added hydrogens; the deposited entries do not",
            );
    let Ok(plan) = PolicySpace::new(AnalysisPolicy::default())
        .vary(dimension)
        .plan()
    else {
        panic!("a one-dimension plan");
    };
    assert!(plan.balanced());
    let [decision] = plan.decisions() else {
        panic!("one decision");
    };
    assert_eq!(decision.class, crate::UncertaintyClass::Interpretive);
    assert!(decision.rationale.contains("modelled"));
    assert!(decision.evidence.contains("PDBbind"));
}

#[test]
fn a_plan_that_dropped_every_run_of_a_level_still_knows_it_is_not_the_whole_product() {
    use molframe_core::contract::SymmetryPolicy;
    let space = PolicySpace::new(AnalysisPolicy::default()).vary(PolicyDimension::symmetry([
        SymmetryPolicy::None,
        SymmetryPolicy::Crystallographic,
    ]));
    let Ok(plan) = space.plan_constrained() else {
        panic!("the default unit survives with no symmetry");
    };
    assert_eq!(plan.cost(), 1);
    assert_eq!(plan.skipped(), 1);
    // The runs alone show one level of one decision, which looks whole.
    assert!(!plan.decompose(|_, _| 0.0).balanced);
}
