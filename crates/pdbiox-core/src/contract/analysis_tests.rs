use super::*;
use crate::diagnostic::{Code, Diagnostic};

#[test]
#[allow(clippy::float_cmp, reason = "reads back a stored value")]
fn a_result_dereferences_to_its_value_so_the_common_path_stays_short() {
    let policy = AnalysisPolicy::default();
    let radius = Analysis::complete(14.2_f64, Coverage::complete(1_960), &policy);
    assert_eq!(*radius, 14.2);
    assert_eq!(radius.status, Status::Complete);
}

#[test]
#[allow(clippy::float_cmp, reason = "an exact ratio of exact counts")]
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
        Coverage::default().fraction(),
        1.0,
        "nothing intended is fully covered"
    );
}

#[test]
fn declining_to_answer_is_a_status_rather_than_an_error() {
    let policy = AnalysisPolicy::default();
    let refused = Analysis::indeterminate(0.0_f64, Coverage::default(), &policy);
    assert_eq!(refused.status, Status::Indeterminate);
    assert!(!refused.status.is_usable());
    assert!(Status::Partial.is_usable());
}

#[test]
fn a_silent_default_that_could_matter_is_the_one_worth_surfacing() {
    let policy = AnalysisPolicy::default();
    let result = Analysis::complete((), Coverage::complete(1), &policy)
        .with_assumption(Assumption::new(
            PolicyField::Hydrogens,
            "ExplicitOnly",
            AssumptionSource::ProfileDefault,
            ImpactEstimate::Moderate,
        ))
        .with_assumption(Assumption::new(
            PolicyField::Altloc,
            "ConformerConsistent",
            AssumptionSource::ProfileDefault,
            ImpactEstimate::None,
        ))
        .with_assumption(Assumption::new(
            PolicyField::Model,
            "First",
            AssumptionSource::Explicit,
            ImpactEstimate::High,
        ));

    let surfaced: Vec<_> = result
        .silent_assumptions()
        .map(|assumption| assumption.field)
        .collect();
    assert_eq!(surfaced, [PolicyField::Hydrogens]);
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
    assert_eq!(*result, 3);
    assert_eq!(result.warnings.len(), 1);
    assert_eq!(result.status, Status::Partial);
}

#[test]
fn mapping_the_value_keeps_everything_that_qualifies_it() {
    let policy = AnalysisPolicy::default();
    let counted = Analysis::partial(vec![1_u32, 2, 3], Coverage::complete(3), &policy)
        .with_warning(Diagnostic::new(Code::W3301))
        .map(|contacts| contacts.len());
    assert_eq!(*counted, 3);
    assert_eq!(counted.status, Status::Partial);
    assert_eq!(counted.warnings.len(), 1);
}
