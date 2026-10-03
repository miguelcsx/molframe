use super::command::sse_name;
use molframe::SecondaryStructure as Ss;

#[test]
fn every_secondary_structure_state_has_a_distinct_stable_name() {
    let names = [
        (Ss::Unknown, "unknown"),
        (Ss::Coil, "coil"),
        (Ss::AlphaHelix, "alpha-helix"),
        (Ss::ThreeTenHelix, "3-10-helix"),
        (Ss::PiHelix, "pi-helix"),
        (Ss::OtherHelix, "other-helix"),
        (Ss::BetaBridge, "beta-bridge"),
        (Ss::Strand, "strand"),
        (Ss::Turn, "turn"),
        (Ss::Bend, "bend"),
    ];
    for (state, name) in names {
        assert_eq!(sse_name(state), name);
    }
}
