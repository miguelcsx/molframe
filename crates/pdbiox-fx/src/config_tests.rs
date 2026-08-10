use super::*;

#[test]
fn complete_toml_specification_builds_native_types() {
    let file = match tempfile::Builder::new().suffix(".toml").tempfile() {
        Ok(file) => file,
        Err(error) => panic!("temporary file failed: {error}"),
    };
    let text = r#"
[components.site]
role = "residue"
component_ids = ["SER"]
required_atoms = ["OG"]

[[constraints]]
kind = "distance"
name = "bond"
first = "site.OG"
second = "site.CB"
target = 1.4
tolerance = 0.2

[profile]
id = "example-1.0"
missing = "indeterminate"

[[profile.rules]]
metric = "bond"
comparison = "between"
values = [1.2, 1.6]
"#;
    if let Err(error) = std::fs::write(file.path(), text) {
        panic!("specification write failed: {error}")
    }
    let specification = match read_evaluation_specification(file.path()) {
        Ok(specification) => specification,
        Err(error) => panic!("specification read failed: {error}"),
    };
    assert_eq!(specification.motif.constraints().len(), 1);
    assert_eq!(specification.profile.id(), "example-1.0");
}
