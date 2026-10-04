use super::*;
use crate::{CategoricalFlip, PolicyDimension, PolicySpace, ScalarError, ScalarMode, SetOverlap};
use molframe_core::contract::{
    AltlocPolicy, AssemblyChoice, Coverage, HydrogenPolicy, Indeterminacy, MissingPolicy,
    SymmetryPolicy,
};
use std::collections::BTreeSet;
use std::convert::Infallible;

fn plan() -> AuditPlan {
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::hydrogens([
            HydrogenPolicy::ExplicitOnly,
            HydrogenPolicy::Exclude,
        ]))
        .vary(PolicyDimension::altloc([
            AltlocPolicy::ConformerConsistent,
            AltlocPolicy::First,
        ]));
    match space.plan() {
        Ok(plan) => plan,
        Err(error) => panic!("plan failed: {error}"),
    }
}

fn close(left: f64, right: f64) -> bool {
    (left - right).abs() < 1e-12
}

fn excluded(policy: &AnalysisPolicy) -> bool {
    matches!(policy.hydrogens, HydrogenPolicy::Exclude)
}

fn first(policy: &AnalysisPolicy) -> bool {
    matches!(policy.altloc, AltlocPolicy::First)
}

#[test]
fn additive_decisions_have_main_effects_and_no_interaction() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let outcome = |policy: &AnalysisPolicy| -> Result<f64, Infallible> {
        Ok(if excluded(policy) { 10.0 } else { 0.0 } + if first(policy) { 1.0 } else { 0.0 })
    };
    let Ok(audit) = audit_outcomes(&plan(), outcome, &metric);
    assert_eq!(audit.metric, "absolute-error");
    let effects = &audit.decomposition.main_effects;
    assert_eq!(effects[0].field, PolicyField::Hydrogens);
    assert!(close(effects[0].mean_change, 10.0));
    assert!(close(effects[1].mean_change, 1.0));
    assert!(close(audit.decomposition.interactions[0].share, 0.0));
    assert!(close(audit.decomposition.higher_order, 0.0));
    // Only the run that takes both first alternatives equals the first run.
    assert!(close(audit.agreement_with_first(), 0.25));
}

#[test]
fn a_pair_that_matters_only_together_is_an_interaction_not_two_small_effects() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let outcome = |policy: &AnalysisPolicy| -> Result<f64, Infallible> {
        Ok(if excluded(policy) && first(policy) {
            10.0
        } else {
            0.0
        })
    };
    let Ok(audit) = audit_outcomes(&plan(), outcome, &metric);
    assert!(close(audit.decomposition.interactions[0].share, 1.0 / 3.0));
    assert!(close(audit.decomposition.main_effects[0].share, 1.0 / 3.0));
}

#[test]
fn the_measured_impact_is_attached_to_the_assumption_it_belongs_to() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let outcome = |policy: &AnalysisPolicy| -> Result<f64, Infallible> {
        Ok(if excluded(policy) { 10.0 } else { 0.0 })
    };
    let Ok(audit) = audit_outcomes(&plan(), outcome, &metric);
    let mut assumptions = vec![
        Assumption::new(
            PolicyField::Hydrogens,
            "explicit-only",
            molframe_core::contract::AssumptionSource::ProfileDefault,
            Impact::Unmeasured,
        ),
        Assumption::new(
            PolicyField::Model,
            "first",
            molframe_core::contract::AssumptionSource::ProfileDefault,
            Impact::Unmeasured,
        ),
    ];
    audit.measure(&mut assumptions);
    assert!(matches!(
        assumptions[0].impact,
        Impact::Measured(measured)
            if close(measured.change, 10.0) && measured.metric == "absolute-error" && measured.comparisons == 2
    ));
    // The model was never varied, so nothing was measured and nothing is claimed.
    assert_eq!(assumptions[1].impact, Impact::Unmeasured);
    assert!(assumptions[0].is_silent_and_material(1.0));
}

#[test]
fn a_conclusion_that_flips_is_counted_not_averaged() {
    let metric = CategoricalFlip::new(|score: &f64| *score >= 0.23);
    let outcome = |policy: &AnalysisPolicy| -> Result<f64, Infallible> {
        Ok(if excluded(policy) { 0.20 } else { 0.64 })
    };
    let Ok(audit) = audit_outcomes(&plan(), outcome, &metric);
    assert!(close(audit.decomposition.main_effects[0].mean_change, 1.0));
    assert!(close(audit.decomposition.main_effects[1].mean_change, 0.0));
    assert!(close(audit.agreement_with_first(), 0.5));
}

fn read_all(policy: &AnalysisPolicy, value: f64, fields: &[PolicyField]) -> Analysis<f64> {
    let mut result = Analysis::complete(value, Coverage::complete(1), policy);
    result.provenance = result.provenance.with_policy_reads(fields);
    result
}

#[test]
fn an_analysis_that_never_read_a_varied_field_refuses_the_audit() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let outcome = |policy: &AnalysisPolicy| -> Result<Analysis<f64>, Infallible> {
        // Reads hydrogens, but not altloc: the second dimension would report stability for free.
        Ok(read_all(
            policy,
            if excluded(policy) { 1.0 } else { 0.0 },
            &[PolicyField::Hydrogens],
        ))
    };
    let refused = audit_analyses(&plan(), outcome, &metric);
    assert!(matches!(
        refused,
        Err(AuditError::Plan(PlanError::NotRead(PolicyField::Altloc)))
    ));
}

#[test]
fn a_later_universe_cannot_borrow_policy_reads_from_the_first() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let mut calls = 0;
    let outcome = |policy: &AnalysisPolicy| -> Result<Analysis<f64>, Infallible> {
        calls += 1;
        let fields = if calls == 1 {
            &[PolicyField::Hydrogens, PolicyField::Altloc][..]
        } else {
            &[PolicyField::Hydrogens][..]
        };
        Ok(read_all(policy, 1.0, fields))
    };
    assert!(matches!(
        audit_analyses(&plan(), outcome, &metric),
        Err(AuditError::Plan(PlanError::NotRead(PolicyField::Altloc)))
    ));
    assert_eq!(calls, 2);
}

#[test]
fn a_later_universe_without_a_read_record_refuses_the_audit() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let mut calls = 0;
    let outcome = |policy: &AnalysisPolicy| -> Result<Analysis<f64>, Infallible> {
        calls += 1;
        Ok(if calls == 1 {
            read_all(policy, 1.0, &[PolicyField::Hydrogens, PolicyField::Altloc])
        } else {
            Analysis::complete(1.0, Coverage::complete(1), policy)
        })
    };
    assert!(matches!(
        audit_analyses(&plan(), outcome, &metric),
        Err(AuditError::Plan(PlanError::NotRead(PolicyField::Hydrogens)))
    ));
    assert_eq!(calls, 2);
}

#[test]
fn universes_with_no_answer_are_a_fraction_and_withhold_the_decomposition() {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let space = PolicySpace::new(AnalysisPolicy::default()).vary(PolicyDimension::missing_atoms([
        MissingPolicy::Report,
        MissingPolicy::Indeterminate,
    ]));
    let Ok(plan) = space.plan() else {
        panic!("a one-dimension plan");
    };
    let outcome = |policy: &AnalysisPolicy| -> Result<Analysis<f64>, Infallible> {
        let mut result = if matches!(policy.missing_atoms, MissingPolicy::Indeterminate) {
            Analysis::indeterminate(
                Indeterminacy::MissingInputs {
                    missing: 2,
                    ambiguous: 0,
                },
                Coverage::default(),
                policy,
            )
        } else {
            Analysis::complete(1.0, Coverage::complete(1), policy)
        };
        result.provenance = result
            .provenance
            .with_policy_reads(&[PolicyField::MissingAtoms]);
        Ok(result)
    };
    let Ok(audit) = audit_analyses(&plan, outcome, &metric) else {
        panic!("a plan over a field the analysis reads");
    };
    assert_eq!(audit.indeterminate, vec![1]);
    assert!(close(audit.indeterminate_fraction(), 0.5));
    assert!(audit.decomposition.is_none());
    assert!(audit.agreement_with_first().is_none());
    assert_eq!(audit.read, vec![PolicyField::MissingAtoms]);
}

#[test]
fn a_space_with_a_contradictory_run_is_refused_before_anything_runs() {
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::assembly([
            AssemblyChoice::AsymmetricUnit,
            AssemblyChoice::Biological("1".into()),
        ]))
        .vary(PolicyDimension::symmetry([
            SymmetryPolicy::None,
            SymmetryPolicy::Crystallographic,
        ]));
    // Run 1 is the unit with crystallographic copies; a biological assembly with them is run 3.
    assert!(matches!(space.plan(), Err(PlanError::Conflict { run: 1 })));
}

#[test]
fn set_overlap_audits_agree_with_the_set_audit() {
    let metric = SetOverlap::new(|items: &BTreeSet<u32>| items.clone());
    let outcome = |policy: &AnalysisPolicy| -> Result<BTreeSet<u32>, Infallible> {
        Ok(if excluded(policy) {
            [1, 2].into()
        } else {
            [1, 2, 3].into()
        })
    };
    let Ok(audit) = audit_outcomes(&plan(), outcome, &metric);
    let Ok(sets) = crate::audit(&plan(), outcome, |items: &BTreeSet<u32>| items.clone());
    assert!(close(
        audit.decomposition.main_effects[0].mean_change,
        sets.dimensions[0].mean_change
    ));
    assert!(close(
        audit.decomposition.main_effects[0].share,
        sets.dimensions[0].main_effect_share
    ));
}

#[test]
fn a_changed_or_omitted_estimand_refuses_before_attribution() {
    for later in [Some("crystal packing interface"), None] {
        let mut calls = 0;
        let result = audit_analyses(
            &plan(),
            |policy| -> Result<Analysis<f64>, Infallible> {
                calls += 1;
                let mut result =
                    read_all(policy, 0.0, &[PolicyField::Hydrogens, PolicyField::Altloc]);
                let estimand = if calls == 1 {
                    Some("biological interface")
                } else {
                    later
                };
                if let Some(text) = estimand {
                    result.provenance = result.provenance.with_estimand(text);
                }
                Ok(result)
            },
            &ScalarError::new(|value: &f64| *value, ScalarMode::Absolute),
        );
        assert!(matches!(
            result,
            Err(AuditError::IncompatibleEstimand { run: 1 })
        ));
        assert_eq!(calls, 2);
    }
}

#[test]
fn an_estimand_cannot_appear_only_after_the_first_universe() {
    let mut calls = 0;
    let result = audit_analyses(
        &plan(),
        |policy| -> Result<Analysis<f64>, Infallible> {
            calls += 1;
            let mut result = read_all(policy, 0.0, &[PolicyField::Hydrogens, PolicyField::Altloc]);
            if calls > 1 {
                result.provenance = result.provenance.with_estimand("ligand contacts");
            }
            Ok(result)
        },
        &ScalarError::new(|value: &f64| *value, ScalarMode::Absolute),
    );
    assert!(matches!(
        result,
        Err(AuditError::IncompatibleEstimand { run: 1 })
    ));
}

#[test]
fn consistent_estimands_keep_the_conclusion_comparable() {
    let result = audit_analyses(
        &plan(),
        |policy| -> Result<Analysis<f64>, Infallible> {
            let mut result = read_all(policy, 0.0, &[PolicyField::Hydrogens, PolicyField::Altloc]);
            result.provenance = result.provenance.with_estimand("ligand contacts");
            Ok(result)
        },
        &ScalarError::new(|value: &f64| *value, ScalarMode::Absolute),
    )
    .unwrap();
    assert!(result.decomposition.is_some());
    assert_eq!(result.runs.len(), 4);
}
