use super::*;

#[test]
fn a_record_carries_the_profile_when_the_policy_is_one_unaltered() {
    let record = Provenance::new(&AnalysisPolicy::default());
    assert_eq!(record.profile, Some(ProfileId::DEFAULT));
    assert_eq!(record.molframe_version, env!("CARGO_PKG_VERSION"));
}

#[test]
fn a_record_of_a_modified_policy_names_no_profile() {
    use super::super::policy::Namespace;
    let policy = AnalysisPolicy::default().with_identifiers(Namespace::Explicit);
    assert_eq!(Provenance::new(&policy).profile, None);
}

#[test]
fn the_clock_does_not_change_the_fingerprint() {
    let policy = AnalysisPolicy::default();
    let without = Provenance::new(&policy);
    let with = Provenance::new(&policy).with_timestamp("2026-08-05T19:00:00Z");
    assert_eq!(without.fingerprint(), with.fingerprint());
}

#[test]
fn a_different_input_changes_the_fingerprint() {
    let policy = AnalysisPolicy::default();
    let first = Provenance::new(&policy).with_input_fingerprint(Fingerprint::of(b"one"));
    let second = Provenance::new(&policy).with_input_fingerprint(Fingerprint::of(b"two"));
    assert_ne!(first.fingerprint(), second.fingerprint());
}

#[test]
fn a_different_reference_dataset_changes_the_fingerprint() {
    let policy = AnalysisPolicy::default();
    let mut first = Provenance::new(&policy);
    first.component_version = Some(DictionaryVersion::new("2026-07-15"));
    let mut second = Provenance::new(&policy);
    second.component_version = Some(DictionaryVersion::new("2026-08-01"));
    assert_ne!(first.fingerprint(), second.fingerprint());
}

#[test]
fn a_source_reports_where_it_came_from() {
    let path = SourceRef::path("entry.cif");
    assert_eq!(path.to_string(), "entry.cif");
    assert!(path.as_path().is_some());
    assert_eq!(SourceRef::Memory.to_string(), "(memory)");
    assert!(SourceRef::None.as_path().is_none());
}

#[test]
fn algorithm_parameters_are_sorted_and_fingerprinted() {
    let policy = AnalysisPolicy::default();
    let Some(epsilon) = ParameterValue::finite_float(2.0) else {
        panic!("finite test value must be accepted");
    };
    let first = Provenance::new(&policy)
        .with_algorithm(AlgorithmId::new("diffusion-map", "1"))
        .with_parameter("epsilon", epsilon.clone())
        .with_parameter("time", ParameterValue::Integer(2));
    let second = Provenance::new(&policy)
        .with_algorithm(AlgorithmId::new("diffusion-map", "1"))
        .with_parameter("time", ParameterValue::Integer(2))
        .with_parameter("epsilon", epsilon);
    assert_eq!(first.fingerprint(), second.fingerprint());
    assert!(first.to_json().contains("parameter.epsilon"));
    let changed = first.with_parameter("time", ParameterValue::Integer(3));
    assert_ne!(changed.fingerprint(), second.fingerprint());
}
