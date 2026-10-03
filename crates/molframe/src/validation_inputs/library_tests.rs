use super::*;
use crate::validation_inputs::ValidationInputError;
use molframe_validate::ReferenceError;

fn write(suffix: &str, text: &str) -> tempfile::NamedTempFile {
    let file = match tempfile::Builder::new().suffix(suffix).tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    if let Err(error) = std::fs::write(file.path(), text) {
        panic!("write failed: {error}")
    }
    file
}

const LIBRARY: &str = "\
[library]
id = \"reference\"
version = \"2026-10\"

[[distribution]]
name = \"c-n\"
edges = [-0.1, 0.0, 0.1]
weights = [1.0, 3.0]

[[distribution]]
name = \"rama\"
x_edges = [-180.0, 0.0, 180.0]
y_edges = [-180.0, 0.0, 180.0]
weights = [1.0, 2.0, 3.0, 4.0]
";

#[test]
fn a_toml_library_builds_the_same_set_as_the_constructors() {
    let file = write(".toml", LIBRARY);
    let read = match read_reference_library(file.path()) {
        Ok(library) => library,
        Err(error) => panic!("library read failed: {error}"),
    };
    let expected = ReferenceLibrary::new(
        "reference",
        "2026-10",
        [
            ReferenceDistribution::histogram("c-n", vec![-0.1, 0.0, 0.1], vec![1.0, 3.0])
                .expect("valid histogram"),
            ReferenceDistribution::grid(
                "rama",
                vec![-180.0, 0.0, 180.0],
                vec![-180.0, 0.0, 180.0],
                vec![1.0, 2.0, 3.0, 4.0],
            )
            .expect("valid grid"),
        ],
    )
    .expect("valid library");
    assert_eq!(read, expected);
}

#[test]
fn a_json_library_is_the_same_document() {
    let file = write(
        ".JSON",
        r#"{"library":{"id":"reference","version":"2026-10"},
            "distribution":[{"name":"c-n","edges":[-0.1,0.0,0.1],"weights":[1.0,3.0]}]}"#,
    );
    let library = match read_reference_library(file.path()) {
        Ok(library) => library,
        Err(error) => panic!("library read failed: {error}"),
    };
    assert_eq!(library.id(), "reference");
    assert_eq!(library.version(), "2026-10");
}

#[test]
fn a_distribution_naming_both_forms_is_ambiguous() {
    let file = write(
        ".toml",
        "[library]\nid='a'\nversion='1'\n[[distribution]]\nname='x'\nedges=[0.0,1.0]\nx_edges=[0.0,1.0]\ny_edges=[0.0,1.0]\nweights=[1.0]\n",
    );
    assert!(matches!(
        read_reference_library(file.path()),
        Err(ValidationInputError::InvalidValue { .. })
    ));
}

#[test]
fn half_a_grid_is_refused() {
    let file = write(
        ".toml",
        "[library]\nid='a'\nversion='1'\n[[distribution]]\nname='x'\nx_edges=[0.0,1.0]\nweights=[1.0]\n",
    );
    assert!(matches!(
        read_reference_library(file.path()),
        Err(ValidationInputError::InvalidValue { .. })
    ));
}

#[test]
fn the_kernels_own_refusals_pass_through() {
    let file = write(
        ".toml",
        "[library]\nid='a'\nversion='1'\n[[distribution]]\nname='x'\nedges=[0.0,1.0]\nweights=[0.0]\n",
    );
    assert!(matches!(
        read_reference_library(file.path()),
        Err(ValidationInputError::Reference(ReferenceError::Weight))
    ));
}

#[test]
fn duplicate_distribution_names_are_refused() {
    let file = write(
        ".toml",
        "[library]\nid='a'\nversion='1'\n[[distribution]]\nname='x'\nedges=[0.0,1.0]\nweights=[1.0]\n[[distribution]]\nname='x'\nedges=[0.0,1.0]\nweights=[1.0]\n",
    );
    assert!(matches!(
        read_reference_library(file.path()),
        Err(ValidationInputError::Reference(ReferenceError::Duplicate))
    ));
}

#[test]
fn unknown_keys_and_unknown_suffixes_are_errors() {
    let unknown = write(
        ".toml",
        "[library]\nid='a'\nversion='1'\nextra=1\n[[distribution]]\nname='x'\nedges=[0.0,1.0]\nweights=[1.0]\n",
    );
    assert!(matches!(
        read_reference_library(unknown.path()),
        Err(ValidationInputError::Toml(_))
    ));
    let suffix = write(".yaml", "");
    assert!(matches!(
        read_reference_library(suffix.path()),
        Err(ValidationInputError::UnsupportedFormat)
    ));
}
