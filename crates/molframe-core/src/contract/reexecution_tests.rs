use super::{ReexecutionEnvironment, ReexecutionError, reexecute_from_provenance};
use crate::contract::{AnalysisPolicy, Fingerprint, Namespace, Provenance};

fn provenance(input: &[u8]) -> Provenance {
    let policy = AnalysisPolicy::default().with_identifiers(Namespace::Explicit);
    Provenance::new(&policy).with_input_fingerprint(Fingerprint::of(input))
}

#[test]
fn replay_receives_exact_bytes_and_recorded_policy() {
    let input = b"data_structure";
    let original = provenance(input);
    let replay = reexecute_from_provenance(
        &original,
        input,
        ReexecutionEnvironment::current(),
        |bytes, policy| (Fingerprint::of(bytes), policy.fingerprint()),
    )
    .unwrap_or_else(|error| panic!("re-execution failed: {error}"));
    assert_eq!(replay.value.0, Fingerprint::of(input));
    assert_eq!(replay.value.1, original.policy_fingerprint);
    assert_eq!(replay.provenance.fingerprint(), original.fingerprint());
}

#[test]
fn replay_refuses_different_input_before_calling_analysis() {
    let original = provenance(b"original");
    let result = reexecute_from_provenance(
        &original,
        b"changed",
        ReexecutionEnvironment::current(),
        |_, _| panic!("callback must not run"),
    );
    assert!(matches!(result, Err(ReexecutionError::InputMismatch)));
}

#[test]
fn replay_refuses_an_unavailable_recorded_dictionary() {
    let input = b"data_structure";
    let mut original = provenance(input);
    original.component_version = Some(crate::contract::DictionaryVersion::new("ccd-2026-08"));
    let result = reexecute_from_provenance(
        &original,
        input,
        ReexecutionEnvironment::current(),
        |_, _| panic!("callback must not run"),
    );
    assert!(matches!(
        result,
        Err(ReexecutionError::ComponentVersionMismatch)
    ));
}
