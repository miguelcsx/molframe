use super::command::sse_name;
use pdbiox::analysis::SseKind;

#[test]
fn every_secondary_structure_state_has_a_stable_name() {
    assert_eq!(sse_name(SseKind::AlphaHelix), "alpha-helix");
    assert_eq!(sse_name(SseKind::Strand), "strand");
    assert_eq!(sse_name(SseKind::Turn), "turn");
    assert_eq!(sse_name(SseKind::Coil), "coil");
}
