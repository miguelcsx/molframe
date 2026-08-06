use super::*;
use crate::span::Position;

fn at(offset: u32) -> ByteSpan {
    ByteSpan::new(Position::new(offset, 1, offset + 1), offset + 1)
}

#[test]
fn a_finding_with_no_message_reports_its_registered_cause() {
    let finding = Diagnostic::new(Code::W3011);
    assert_eq!(
        finding.message(),
        "residue boundary was ambiguous and fell back to file order"
    );
    assert!(finding.remedy().contains("residue identifiers"));
}

#[test]
fn a_finding_carries_the_severity_its_code_is_registered_with() {
    assert_eq!(Diagnostic::new(Code::E1103).severity(), Severity::Breaking);
    assert_eq!(Diagnostic::new(Code::W2001).severity(), Severity::Info);
}

#[test]
fn findings_are_returned_worst_first_then_by_position_then_by_code() {
    let mut sink = Diagnostics::new();
    sink.push(Diagnostic::new(Code::W2001).at(at(50)));
    sink.push(Diagnostic::new(Code::E1103).at(at(90)));
    sink.push(Diagnostic::new(Code::W3011).at(at(10)));
    sink.push(Diagnostic::new(Code::W2002).at(at(50)));

    let ordered: Vec<_> = sink.finish().iter().map(Diagnostic::code).collect();
    assert_eq!(
        ordered,
        [Code::E1103, Code::W3011, Code::W2001, Code::W2002]
    );
}

#[test]
fn findings_without_a_position_sort_after_findings_with_one() {
    let mut sink = Diagnostics::new();
    sink.push(Diagnostic::new(Code::W2001));
    sink.push(Diagnostic::new(Code::W2002).at(at(900)));

    let ordered: Vec<_> = sink.finish().iter().map(Diagnostic::code).collect();
    assert_eq!(ordered, [Code::W2002, Code::W2001]);
}

#[test]
fn strictness_decides_whether_an_accumulated_finding_counts_as_failure() {
    let mut sink = Diagnostics::new();
    sink.push(Diagnostic::new(Code::E2003));
    assert!(sink.has_error(Strictness::Strict));
    assert!(!sink.has_error(Strictness::Medium));
    assert!(!sink.has_error(Strictness::Loose));
}

#[test]
fn the_worst_severity_is_tracked_without_rescanning_the_list() {
    let mut sink = Diagnostics::new();
    assert_eq!(sink.worst(), None);
    sink.push(Diagnostic::new(Code::W2001));
    assert_eq!(sink.worst(), Some(Severity::Info));
    sink.push(Diagnostic::new(Code::E1103));
    sink.push(Diagnostic::new(Code::W2002));
    assert_eq!(sink.worst(), Some(Severity::Breaking));
}

#[test]
fn the_same_input_produces_the_same_list_in_the_same_order() {
    let build = || {
        let mut sink = Diagnostics::new();
        for (code, offset) in [
            (Code::W2002, 30),
            (Code::E1101, 5),
            (Code::W3011, 30),
            (Code::W2001, 1),
        ] {
            sink.push(Diagnostic::new(code).at(at(offset)));
        }
        sink.finish()
            .iter()
            .map(Diagnostic::code)
            .collect::<Vec<_>>()
    };
    assert_eq!(build(), build());
}
