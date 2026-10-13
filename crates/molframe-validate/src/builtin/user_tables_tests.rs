use super::*;
use crate::builtin_validation::{ReferenceSet, bundled_reference};

fn scratch(name: &str, bytes: &[u8]) -> Result<std::path::PathBuf, std::io::Error> {
    let dir = std::env::temp_dir().join(format!("molframe-validate-{}", std::process::id()));
    std::fs::create_dir_all(&dir)?;
    let path = dir.join(name);
    std::fs::write(&path, bytes)?;
    Ok(path)
}

#[test]
fn a_user_geometry_file_replaces_the_default_and_is_reported()
-> Result<(), Box<dyn std::error::Error>> {
    let path = scratch("geometry.tsv", b"bond\tALA\tN-CA\t1.5\t0.02\tX\n")?;
    let overrides = ReferenceOverrides::from_paths(None, None, Some(&path), None)?;
    let set: ReferenceSet = bundled_reference()?.clone().overlaid(overrides);
    let info = set.geometry.info();
    assert_eq!(info.name, "user-geometry");
    assert_eq!(info.version, "user-supplied");
    assert!(info.source.ends_with("geometry.tsv"));
    assert_eq!(set.geometry.len(), 1);
    assert_eq!(set.rotamers.info().licence, "CC-BY-4.0");
    Ok(())
}

#[test]
fn a_user_rama_file_round_trips_through_the_decoder() -> Result<(), Box<dyn std::error::Error>> {
    let bytes = std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/data/rama.bin"))?;
    let path = scratch("rama.bin", &bytes)?;
    let overrides = ReferenceOverrides::from_paths(Some(&path), None, None, None)?;
    let part = overrides.rama.ok_or("missing rama")?;
    let default = &bundled_reference()?.rama;
    let category = crate::builtin_validation::RamaCategory::General;
    let decoded = part.density(category, -63.0, -43.0);
    assert!((decoded - default.density(category, -63.0, -43.0)).abs() < 1e-6);
    Ok(())
}

#[test]
fn bad_and_missing_files_are_errors() -> Result<(), Box<dyn std::error::Error>> {
    let bad = scratch("bad.bin", b"nope")?;
    assert!(matches!(
        ReferenceOverrides::from_paths(Some(&bad), None, None, None),
        Err(BuiltinError::Table(_))
    ));
    let empty = scratch("empty.tsv", b"# nothing\n")?;
    assert!(matches!(
        ReferenceOverrides::from_paths(None, None, Some(&empty), None),
        Err(BuiltinError::Load { .. })
    ));
    let missing = std::path::Path::new("/nonexistent/molframe/chi.tsv");
    assert!(matches!(
        ReferenceOverrides::from_paths(None, None, None, Some(missing)),
        Err(BuiltinError::Load { .. })
    ));
    Ok(())
}
