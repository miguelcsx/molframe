use super::*;

const GOOD: &str = r#"
[[entry]]
name = "lab-radii"
version = "2"
url = "https://example.org/radii.tsv"
licence = "CC0-1.0"
sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
extension = "radii"
converter = "tsv"
append = ["extra.rows"]
"#;

fn load(text: &str) -> Result<Manifest, DataError> {
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join("manifest.toml");
    std::fs::write(&path, text).expect("write");
    Manifest::load(&path)
}

#[test]
fn a_valid_manifest_exposes_its_entries_and_cache_names() {
    let manifest = load(GOOD).expect("loads");
    let entry = manifest.entry("lab-radii").expect("entry");
    assert_eq!(entry.file_name(), "lab-radii.radii");
    assert_eq!(entry.append.len(), 1);
    assert!(matches!(
        manifest.entry("other"),
        Err(DataError::NoEntry { .. })
    ));
}

#[test]
fn schema_violations_are_refused_naming_the_entry() {
    for (from, to) in [
        ("https://example.org", "http://example.org"),
        ("000000000000", "zz0000000000"),
        ("name = \"lab-radii\"", "name = \"a/b\""),
        ("licence = \"CC0-1.0\"", "licence = \"\""),
        (
            "extension = \"radii\"",
            "extension = \"radii\"\nsurprise = 1",
        ),
    ] {
        let error = load(&GOOD.replace(from, to)).expect_err("refused");
        assert!(
            matches!(error, DataError::Manifest { .. }),
            "{from}: {error}"
        );
    }
    let twice = format!("{GOOD}{GOOD}");
    assert!(
        load(&twice)
            .expect_err("duplicate")
            .to_string()
            .contains("twice")
    );
}
