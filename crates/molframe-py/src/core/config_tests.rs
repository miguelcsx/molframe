use super::PyPolicyOverrides;

#[test]
fn policy_overrides_round_trip_into_the_native_policy_parser() {
    let overrides = PyPolicyOverrides {
        assembly: Some("asymmetric-unit".to_owned()),
        model: Some("first".to_owned()),
        altloc: None,
        identifiers: Some("auth".to_owned()),
        missing_atoms: Some("report".to_owned()),
        hydrogens: None,
        atom_equivalence: None,
        symmetry: None,
        alignment: None,
        precision: None,
        periodic: None,
        vdw_radii: None,
        contact_def: None,
        float_tolerance_relative: None,
        float_tolerance_absolute: None,
    };
    let policy = overrides
        .into_inner()
        .apply_to(molframe::AnalysisPolicy::default())
        .expect("known policy vocabulary");
    assert_eq!(policy.identifiers, molframe::Namespace::Auth);
}
