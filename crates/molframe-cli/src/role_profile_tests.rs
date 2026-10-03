use super::read;
use molframe::chemistry::{ComponentKind, PolymerAtomRole};

fn profile_file(suffix: &str, contents: &str) -> tempfile::NamedTempFile {
    let file = tempfile::Builder::new()
        .suffix(suffix)
        .tempfile()
        .expect("temporary profile");
    std::fs::write(file.path(), contents).expect("write profile");
    file
}

#[test]
fn json_and_toml_preserve_exact_profile_names_constraints_and_combined_roles() {
    let json = profile_file(
        ".JSON",
        r#"{"profile_id":"caller-v2","rules":[{"atom_name":"Q","role":3,"component_id":"GLY","component_kind":1}]}"#,
    );
    let toml = profile_file(
        ".TOML",
        "profile_id = 'caller-v2'\n[[rules]]\natom_name = 'Q'\nrole = 3\ncomponent_id = 'GLY'\ncomponent_kind = 1\n",
    );
    let profile = read(json.path()).expect("valid JSON");
    assert_eq!(profile, read(toml.path()).expect("valid TOML"));
    assert_eq!(profile.id.as_ref(), "caller-v2");
    assert_eq!(profile.rules[0].atom_name.as_ref(), "Q");
    assert_eq!(profile.rules[0].component_id.as_deref(), Some("GLY"));
    assert_eq!(
        profile.rules[0].component_kind,
        Some(ComponentKind::AminoAcid)
    );
    assert_eq!(
        profile.rules[0].role,
        PolymerAtomRole::PROTEIN_NITROGEN.union(PolymerAtomRole::PROTEIN_ALPHA_CARBON)
    );
}

#[test]
fn invalid_codes_unknown_fields_and_missing_required_fields_are_not_defaulted() {
    for contents in [
        r#"{"profile_id":"v1","rules":[{"atom_name":"N","role":-1,"component_kind":1}]}"#,
        r#"{"profile_id":"v1","rules":[{"atom_name":"N","role":262144,"component_kind":1}]}"#,
        r#"{"profile_id":"v1","rules":[{"atom_name":"N","role":1,"component_kind":8}]}"#,
        r#"{"profile_id":"v1","rules":[{"atom_name":"N","role":1,"component_kind":1,"typo":true}]}"#,
        r#"{"profile_id":"v1","rules":[],"typo":true}"#,
        r#"{"rules":[]}"#,
        r#"{"profile_id":"v1"}"#,
        r#"{"profile_id":"v1","rules":[{"atom_name":"N","component_kind":1}]}"#,
    ] {
        let file = profile_file(".json", contents);
        assert!(read(file.path()).is_err(), "accepted {contents}");
    }
}

#[test]
fn missing_files_and_unsupported_suffixes_are_explicit_errors() {
    let file = profile_file(".yaml", "profile_id: v1\nrules: []");
    assert!(read(file.path()).is_err());
    let directory = tempfile::tempdir().expect("temporary directory");
    assert!(read(&directory.path().join("absent.json")).is_err());
}
