use super::*;
use crate::contract::{MissingPolicy, MissingPolicyError, resolve_missing};
use crate::diagnostic::{Code, Diagnostic};

#[test]
fn the_answer_is_read_through_the_outcome_so_asking_whether_it_exists_is_the_way_in() {
    let policy = AnalysisPolicy::default();
    let radius = Analysis::complete(14.2_f64, Coverage::complete(1_960), &policy);
    assert_eq!(
        radius.value().map(|value| value.to_bits()),
        Some(14.2_f64.to_bits())
    );
    assert_eq!(radius.status(), Status::Complete);
    assert_eq!(radius.quality(), Some(Quality::Complete));
}

#[test]
fn coverage_reports_the_fraction_of_intended_atoms_that_were_used() {
    let partial = Coverage {
        intended: 512,
        used: 481,
        missing: 31,
        ambiguous: 0,
    };
    assert!((partial.fraction() - 0.939_453).abs() < 1e-5);
    assert!(!partial.is_complete());
    assert!(Coverage::complete(10).is_complete());
    assert_eq!(
        Coverage::default().fraction().to_bits(),
        1.0_f64.to_bits(),
        "nothing intended is fully covered"
    );
}

#[test]
fn declining_to_answer_has_no_value_to_misread() {
    let policy = AnalysisPolicy::default();
    let refused = Analysis::<f64>::indeterminate(
        Indeterminacy::MissingInputs {
            missing: 31,
            ambiguous: 0,
        },
        Coverage::default(),
        &policy,
    );
    assert_eq!(refused.status(), Status::Indeterminate);
    assert!(!refused.status().is_usable());
    assert!(Status::Partial.is_usable());
    assert_eq!(refused.value(), None);
    assert_eq!(refused.quality(), None);
    assert_eq!(
        refused.indeterminacy().map(ToString::to_string),
        Some("the policy rejects 31 missing and 0 ambiguous inputs".to_owned())
    );
    assert!(matches!(
        refused.into_result(),
        Err(Indeterminacy::MissingInputs { missing: 31, .. })
    ));
}

#[test]
fn degrading_a_refusal_changes_nothing_and_degrading_an_answer_keeps_the_worst() {
    let policy = AnalysisPolicy::default();
    let mut answer = Analysis::complete(1_u32, Coverage::complete(1), &policy);
    answer.degrade(Quality::Partial);
    answer.degrade(Quality::Complete);
    assert_eq!(answer.status(), Status::Partial);
    answer.degrade(Quality::Ambiguous);
    assert_eq!(answer.status(), Status::Ambiguous);

    let mut refused = Analysis::<u32>::indeterminate(
        Indeterminacy::Other("no defensible answer".into()),
        Coverage::default(),
        &policy,
    );
    refused.degrade(Quality::Ambiguous);
    assert_eq!(refused.status(), Status::Indeterminate);
}

#[test]
fn an_impact_is_a_measurement_or_it_is_unknown() {
    let policy = AnalysisPolicy::default();
    let measured = |change| {
        Impact::Measured(MeasuredImpact {
            metric: "jaccard-distance",
            change,
            comparisons: 6,
        })
    };
    let result = Analysis::complete((), Coverage::complete(1), &policy)
        .with_assumption(Assumption::new(
            PolicyField::Hydrogens,
            "ExplicitOnly",
            AssumptionSource::ProfileDefault,
            measured(0.31),
        ))
        .with_assumption(Assumption::new(
            PolicyField::Altloc,
            "ConformerConsistent",
            AssumptionSource::ProfileDefault,
            measured(0.01),
        ))
        .with_assumption(Assumption::new(
            PolicyField::VdwRadii,
            "Bondi",
            AssumptionSource::ProfileDefault,
            Impact::Unmeasured,
        ))
        .with_assumption(Assumption::new(
            PolicyField::Model,
            "First",
            AssumptionSource::Explicit,
            measured(0.9),
        ));

    // Only a measured, silent change at or above the threshold is worth surfacing; an
    // unmeasured default is unknown, not harmless, and not claimed to be material.
    let surfaced: Vec<_> = result
        .silent_assumptions(0.1)
        .map(|assumption| assumption.field)
        .collect();
    assert_eq!(surfaced, [PolicyField::Hydrogens]);
    assert_eq!(Impact::default(), Impact::Unmeasured);
}

#[test]
fn warnings_travel_with_the_result_rather_than_replacing_it() {
    let policy = AnalysisPolicy::default();
    let result = Analysis::partial(
        3_u32,
        Coverage {
            intended: 4,
            used: 3,
            missing: 1,
            ambiguous: 0,
        },
        &policy,
    )
    .with_warning(Diagnostic::new(Code::W3301));
    assert_eq!(result.value(), Some(&3));
    assert_eq!(result.warnings.len(), 1);
    assert_eq!(result.status(), Status::Partial);
}

#[test]
fn mapping_the_value_keeps_everything_that_qualifies_it() {
    let policy = AnalysisPolicy::default();
    let counted = Analysis::partial(vec![1_u32, 2, 3], Coverage::complete(3), &policy)
        .with_warning(Diagnostic::new(Code::W3301))
        .map(|contacts| contacts.len());
    assert_eq!(counted.value(), Some(&3));
    assert_eq!(counted.status(), Status::Partial);
    assert_eq!(counted.warnings.len(), 1);
}

#[test]
fn quality_combines_to_the_weaker_and_names_its_status() {
    assert_eq!(Quality::Complete.worst(Quality::Partial), Quality::Partial);
    assert_eq!(
        Quality::Partial.worst(Quality::Ambiguous),
        Quality::Ambiguous
    );
    assert_eq!(
        Quality::Complete.worst(Quality::Complete),
        Quality::Complete
    );
    assert_eq!(Status::from(Quality::Ambiguous), Status::Ambiguous);
}

#[test]
fn an_indeterminate_frame_explains_the_series() {
    let reason = Indeterminacy::Frame {
        frame: 3,
        reason: Box::new(Indeterminacy::UnresolvedConformations),
    };
    assert_eq!(
        reason.to_string(),
        "frame 3: the alternate-conformation rule selected no consistent atoms"
    );
}

#[test]
fn mapping_a_refusal_leaves_the_refusal() {
    let outcome: Outcome<u32> = Outcome::Indeterminate(Indeterminacy::Other("no".into()));
    assert!(!outcome.map(|value| value + 1).is_determinate());
    assert_eq!(
        Outcome::Determinate(2_u32).map(|value| value + 1).value(),
        Some(&3)
    );
}

#[test]
fn the_missing_policy_decides_what_an_incomplete_value_becomes() {
    let coverage = Coverage {
        intended: 10,
        used: 8,
        missing: 2,
        ambiguous: 0,
    };
    let decide = |policy| resolve_missing(7_u32, Quality::Complete, coverage, policy);
    assert!(matches!(
        decide(MissingPolicy::Ignore),
        Ok((Outcome::Determinate(7), Quality::Complete))
    ));
    assert!(matches!(
        decide(MissingPolicy::Report),
        Ok((Outcome::Determinate(7), Quality::Partial))
    ));
    assert!(matches!(
        decide(MissingPolicy::Indeterminate),
        Ok((
            Outcome::Indeterminate(Indeterminacy::MissingInputs { missing: 2, .. }),
            _
        ))
    ));
    assert_eq!(
        decide(MissingPolicy::Fail).err(),
        Some(MissingPolicyError::Fail {
            missing: 2,
            ambiguous: 0
        })
    );
    // Nothing missing: every policy keeps the value.
    let whole = Coverage::complete(10);
    for policy in [MissingPolicy::Fail, MissingPolicy::Indeterminate] {
        assert!(matches!(
            resolve_missing(7_u32, Quality::Complete, whole, policy),
            Ok((Outcome::Determinate(7), Quality::Complete))
        ));
    }
}
