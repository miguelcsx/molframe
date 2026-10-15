use super::*;

/// The tables the feature compiles in, as `data/` file names.
const EMBEDDED_FILES: [&str; 4] = [
    "macromolecular-reach.reach",
    "pyykko-covalent.covalent",
    "bondi.radii",
    "macromolecular-metals.metals",
];

fn read(path: &std::path::Path) -> Vec<u8> {
    match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(error) => panic!("{}: {error}", path.display()),
    }
}

fn root_data() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data")
}

fn beside_data() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data")
}

#[test]
fn every_embedded_table_beside_the_crate_is_the_repository_copy() {
    for file in EMBEDDED_FILES {
        let beside = beside_data().join(file);
        assert!(beside.is_file(), "{file} is not copied beside the crate");
        assert_eq!(
            read(&beside),
            read(&root_data().join(file)),
            "{file} drifted"
        );
    }
}

#[cfg(feature = "embedded-tables")]
#[test]
fn every_compiled_entry_is_the_repository_file() {
    for (name, extension, bytes) in ENTRIES {
        let file = format!("{name}.{extension}");
        assert!(
            EMBEDDED_FILES.contains(&file.as_str()),
            "{file} is unlisted"
        );
        assert_eq!(bytes, read(&root_data().join(&file)).as_slice(), "{file}");
    }
    assert_eq!(ENTRIES.len(), EMBEDDED_FILES.len());
}

/// Decodes the entry the bundled profile names for `T`, checking the name.
#[cfg(feature = "embedded-tables")]
fn bundled_table<T: crate::TableKind + Send + Sync + 'static>(
    expected: &str,
) -> Arc<crate::Provided<T>> {
    let Some(tables) = bundled() else {
        panic!("the feature compiles the layer in");
    };
    let Some(name) = tables.profile_entry(T::KIND) else {
        panic!("the profile names a {} table", T::KIND);
    };
    assert_eq!(&*name, expected);
    let Some(entry) = tables.entry(&name, T::EXTENSION) else {
        panic!("the {} entry is held", T::KIND);
    };
    match entry.provided::<T>() {
        Ok(provided) => provided,
        Err(error) => panic!("the {} table decodes: {error}", T::KIND),
    }
}

#[cfg(feature = "embedded-tables")]
#[test]
fn the_bundled_profile_names_a_reach_table_that_decodes() {
    let provided = bundled_table::<crate::BondReach>("macromolecular-reach");
    assert_eq!(&provided.table, crate::test_data::reach());
}

#[cfg(feature = "embedded-tables")]
#[test]
fn the_bundled_profile_names_the_tables_the_default_flows_read() {
    let covalent = bundled_table::<crate::CovalentRadii>("pyykko-covalent");
    let shipped = crate::test_data::shipped::<crate::CovalentRadii>("pyykko-covalent.covalent");
    assert_eq!(covalent.table, shipped.table);
    let vdw = bundled_table::<crate::VdwRadii>("bondi");
    assert_eq!(
        vdw.table,
        crate::test_data::shipped::<crate::VdwRadii>("bondi.radii").table
    );
    let sites = bundled_table::<crate::MetalSites>("macromolecular-metals");
    assert_eq!(
        sites.table,
        crate::test_data::shipped::<crate::MetalSites>("macromolecular-metals.metals").table
    );
}

#[cfg(feature = "embedded-tables")]
#[test]
fn the_bundled_layer_is_built_once() {
    let (Some(first), Some(second)) = (bundled(), bundled()) else {
        panic!("the feature compiles the layer in");
    };
    assert!(Arc::ptr_eq(&first, &second));
}

#[cfg(not(feature = "embedded-tables"))]
#[test]
fn without_the_feature_nothing_is_compiled_in() {
    assert!(bundled().is_none());
}

#[cfg(feature = "embedded-tables")]
#[path = "table_bundled_store_tests.rs"]
mod store;
