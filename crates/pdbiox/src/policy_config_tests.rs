use super::*;

#[test]
fn strict_toml_document_applies_named_overrides() {
    let file = match tempfile::Builder::new().suffix(".toml").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary policy failed: {error}"),
    };
    if let Err(error) = std::fs::write(
        file.path(),
        "[policy]\nassembly='biological:1'\nmodel='index:2'\naltloc='label:A'\nidentifiers='explicit'\nmissing_atoms='fail'\nhydrogens='exclude'\natom_equivalence='explicit'\nsymmetry='crystallographic'\nalignment='global'\nprecision='f32'\nperiodic='minimum-image'\nvdw_radii='alvarez'\ncontact_def='surface:1.4'\nfloat_tolerance_relative=0.000001\nfloat_tolerance_absolute=0.000000001\n",
    ) {
        panic!("policy write failed: {error}")
    }
    let policy = match read_policy(file.path()) {
        Ok(policy) => policy,
        Err(error) => panic!("policy read failed: {error}"),
    };
    assert_eq!(policy.model, ModelChoice::Index(2));
    assert_eq!(policy.identifiers, Namespace::Explicit);
    assert_eq!(policy.vdw_radii, RadiiSet::Alvarez);
    let ContactDefinition::SurfaceBased { probe } = policy.contact_def else {
        panic!("surface contact policy absent")
    };
    assert!((probe - 1.4).abs() < f32::EPSILON);
}

#[test]
fn unknown_fields_are_rejected() {
    let file = match tempfile::Builder::new().suffix(".json").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary policy failed: {error}"),
    };
    if let Err(error) = std::fs::write(file.path(), r#"{"policy":{"magic":true}}"#) {
        panic!("policy write failed: {error}")
    }
    assert!(matches!(
        read_policy(file.path()),
        Err(PolicyConfigError::Json(_))
    ));
}

#[test]
fn application_configuration_accepts_output_and_chemistry_sections() {
    let file = match tempfile::Builder::new().suffix(".toml").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary configuration failed: {error}"),
    };
    if let Err(error) = std::fs::write(
        file.path(),
        "[policy]\nmodel='first'\n[output]\nformat='jsonl'\n[chem]\nccd_cache='/data/ccd'\nccd_version='2026-08-01'\n",
    ) {
        panic!("configuration write failed: {error}")
    }
    let configuration = match read_configuration(file.path()) {
        Ok(configuration) => configuration,
        Err(error) => panic!("configuration read failed: {error}"),
    };
    assert_eq!(configuration.output.format.as_deref(), Some("jsonl"));
    assert_eq!(
        configuration.chem.ccd_cache.as_deref(),
        Some(std::path::Path::new("/data/ccd"))
    );
    assert_eq!(
        configuration.chem.ccd_version.as_deref(),
        Some("2026-08-01")
    );
}
