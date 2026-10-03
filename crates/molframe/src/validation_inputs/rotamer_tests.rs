use super::*;
use crate::validation_inputs::ValidationInputError;
use molframe_validate::RotamerError;

fn write(text: &str) -> tempfile::NamedTempFile {
    let file = match tempfile::Builder::new().suffix(".toml").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = std::fs::write(file.path(), text) {
        panic!("write failed: {error}")
    }
    file
}

const PROFILE: &str = "\
[profile]
id = \"chi-wells\"
version = \"2026-10\"

[[definition]]
component_id = \"SER\"
chi_index = 1
atoms = [\"N\", \"CA\", \"CB\", \"OG\"]
distribution = \"ser-chi1\"
";

#[test]
fn a_profile_round_trips_into_the_kernels_type() {
    let file = write(PROFILE);
    let read = match read_rotamer_profile(file.path()) {
        Ok(profile) => profile,
        Err(error) => panic!("profile read failed: {error}"),
    };
    let expected = RotamerProfile::new(
        "chi-wells",
        "2026-10",
        [RotamerDefinition {
            component_id: "SER".into(),
            chi_index: 1,
            atoms: ["N".into(), "CA".into(), "CB".into(), "OG".into()],
            distribution: "ser-chi1".into(),
        }],
    )
    .expect("valid profile");
    assert_eq!(read, expected);
}

#[test]
fn a_torsion_with_three_atoms_is_a_schema_error() {
    let file =
        write(&PROFILE.replace("[\"N\", \"CA\", \"CB\", \"OG\"]", "[\"N\", \"CA\", \"CB\"]"));
    assert!(matches!(
        read_rotamer_profile(file.path()),
        Err(ValidationInputError::Toml(_))
    ));
}

#[test]
fn a_repeated_torsion_is_refused_by_the_kernel() {
    let doubled = format!(
        "{PROFILE}\n{}",
        &PROFILE[PROFILE.find("[[definition]]").unwrap_or(0)..]
    );
    let file = write(&doubled);
    assert!(matches!(
        read_rotamer_profile(file.path()),
        Err(ValidationInputError::Rotamer(RotamerError::InvalidProfile))
    ));
}

#[test]
fn a_profile_without_definitions_is_refused() {
    let file = write("definition = []\n[profile]\nid = 'a'\nversion = '1'\n");
    assert!(matches!(
        read_rotamer_profile(file.path()),
        Err(ValidationInputError::Rotamer(RotamerError::InvalidProfile))
    ));
}
