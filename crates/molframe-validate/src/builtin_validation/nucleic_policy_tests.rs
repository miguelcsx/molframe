use super::*;
use crate::builtin_validation::fixtures::{NUCLEIC, identity};

fn decode(text: &str) -> Result<Provided<NucleicPolicy>, TableTextError> {
    NucleicPolicy::decode(text.as_bytes(), identity("fallback"))
}

#[test]
fn windows_wrap_around_zero() {
    assert!(within(350.0, [-18.0, 54.0]));
    assert!(within(10.0, [-18.0, 54.0]));
    assert!(!within(100.0, [-18.0, 54.0]));
}

#[test]
fn a_policy_table_gives_every_window() -> Result<(), TableTextError> {
    let provided = decode(NUCLEIC)?;
    let policy = provided.table;
    assert_eq!(provided.identity.name, "synthetic-nucleic");
    assert!(
        policy
            .c3_endo_phase
            .iter()
            .zip([-18.0, 54.0])
            .all(|(a, b)| (a - b).abs() < 1e-9)
    );
    assert_eq!(policy.gamma.len(), 3);
    assert_eq!(policy.pucker_class(10.0), PUCKER_C3_ENDO);
    assert_eq!(policy.pucker_class(160.0), PUCKER_C2_ENDO);
    assert_eq!(policy.pucker_class(90.0), PUCKER_OTHER);
    assert!(policy.delta_agrees(PUCKER_C3_ENDO, 80.0));
    assert!(!policy.delta_agrees(PUCKER_C3_ENDO, 150.0));
    Ok(())
}

#[test]
fn an_incomplete_or_malformed_policy_is_refused() {
    let head = "name: x\nversion: 1\nsource: s\nlicence: CC0-1.0\n";
    let windows =
        "c3_endo_phase: 0 1\nc2_endo_phase: 0 1\nc3_endo_delta: 0 1\nc2_endo_delta: 0 1\n";
    assert!(decode(&format!("{head}{windows}")).is_err(), "no gamma row");
    assert!(decode(&format!("{head}30 90\n")).is_err(), "no windows");
    assert!(decode(&format!("{head}{windows}30\n")).is_err());
    assert!(decode(&format!("{head}{windows}c3_endo_phase: 0 1\n30 90\n")).is_err());
    assert!(decode(&format!("{head}{windows}30 x\n")).is_err());
}
