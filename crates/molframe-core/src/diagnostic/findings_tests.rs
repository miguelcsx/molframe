//! A finished set of findings, and the error it makes.

use super::*;
use crate::diagnostic::Diagnostics;

#[test]
fn a_set_reports_the_code_that_stopped_it_rather_than_the_first_one_pushed() {
    let mut diagnostics = Diagnostics::new();
    diagnostics.push(Diagnostic::new(Code::W1002));
    diagnostics.push(Diagnostic::new(Code::E1001));

    // `code()` reads the finished order, so the error leads even though the
    // warning was raised first.
    let findings = Findings::from(diagnostics.finish());
    assert_eq!(findings.code(), Some(Code::E1001));
}

#[test]
fn a_set_with_nothing_in_it_stopped_nothing() {
    let findings = Findings::from(Vec::new());
    assert!(findings.is_empty());
    assert_eq!(findings.code(), None);
}

#[test]
fn a_set_prints_its_first_finding_and_counts_the_rest() {
    let findings = Findings::from(vec![
        Diagnostic::new(Code::E1001),
        Diagnostic::new(Code::W1002),
        Diagnostic::new(Code::W1002),
    ]);
    assert_eq!(
        findings.to_string(),
        format!("{} (and 2 more)", Diagnostic::new(Code::E1001))
    );
}

#[test]
fn a_lone_finding_prints_as_itself() {
    let findings = Findings::from(Diagnostic::new(Code::E1001));
    assert_eq!(
        findings.to_string(),
        Diagnostic::new(Code::E1001).to_string()
    );
}

#[test]
fn a_set_iterates_by_reference_and_derefs_to_its_findings() {
    let findings = Findings::from(vec![
        Diagnostic::new(Code::E1001),
        Diagnostic::new(Code::W1002),
    ]);

    let mut codes = Vec::new();
    for finding in &findings {
        codes.push(finding.code());
    }
    assert_eq!(codes, vec![Code::E1001, Code::W1002]);

    // The slice vocabulary, through `Deref`, without re-declaring any of it.
    assert_eq!(findings.first().map(Diagnostic::code), Some(Code::E1001));
    assert_eq!(findings.as_slice().len(), 2);
}

#[test]
fn a_set_survives_moving_in_and_out_of_a_vector() {
    let findings = vec![Diagnostic::new(Code::E1001), Diagnostic::new(Code::W1002)];
    assert_eq!(Findings::from(findings.clone()).into_vec(), findings);
}

#[test]
fn a_set_is_an_error_a_question_mark_carries() {
    fn stopped() -> Result<(), Box<dyn std::error::Error>> {
        Err(Findings::from(Diagnostic::new(Code::E1001)))?;
        Ok(())
    }

    let Err(error) = stopped() else {
        panic!("a set of findings is an error, so `?` should have returned it")
    };
    assert_eq!(error.to_string(), Diagnostic::new(Code::E1001).to_string());
}
