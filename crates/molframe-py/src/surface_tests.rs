//! Gates on the Python surface, as ordinary tests so `cargo test` enforces them.
//!
//! The extension module is built in-process and compared with what the package
//! promises: the hand-written stubs, the operation catalog, the GIL discipline
//! of every `#[pyfunction]`, the size of the Python sources and the crate's
//! dependency list. Each check names what it found, never just that it failed.
//!
//! Every check needs the whole surface registered, so they run with `full`.

use pyo3::prelude::*;
use pyo3::types::{PyCFunction, PyModule, PyType};
use std::collections::BTreeSet;
use std::path::Path;

#[path = "surface_tests/sources_tests.rs"]
mod sources;
#[path = "surface_tests/stubs_tests.rs"]
mod stubs;
use stubs::{
    CONTRACT_DUNDERS, PYTHON_DEFINED, UNCATALOGUED_NAMESPACES, dunder_all, parse_stub, python_root,
    read, reexports,
};

fn with_native<T>(check: impl FnOnce(&Bound<'_, PyModule>) -> T) -> T {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "_native") {
            Ok(module) => module,
            Err(error) => panic!("the test module should be created: {error}"),
        };
        if let Err(error) = crate::module::native(&module) {
            panic!("the extension module should register: {error}");
        }
        check(&module)
    })
}

fn public_attributes<'py>(object: &Bound<'py, PyAny>) -> Vec<(String, Bound<'py, PyAny>)> {
    let Ok(listing) = object.dir() else {
        panic!("dir() should work on a registered object")
    };
    listing
        .iter()
        .filter_map(|entry| entry.extract::<String>().ok())
        .filter(|name| !name.starts_with('_'))
        .filter_map(|name| {
            object
                .getattr(name.as_str())
                .ok()
                .map(|value| (name, value))
        })
        .collect()
}

/// The registered namespaces: every public attribute of the module that is a module.
fn namespaces<'py>(module: &Bound<'py, PyModule>) -> Vec<(String, Bound<'py, PyAny>)> {
    public_attributes(module.as_any())
        .into_iter()
        .filter(|(_, value)| value.is_instance_of::<PyModule>())
        .collect()
}

fn non_module_names(object: &Bound<'_, PyAny>) -> BTreeSet<String> {
    public_attributes(object)
        .into_iter()
        .filter(|(_, value)| !value.is_instance_of::<PyModule>())
        .map(|(name, _)| name)
        .collect()
}

fn native_members(class: &Bound<'_, PyAny>) -> BTreeSet<String> {
    let Ok(listing) = class.dir() else {
        panic!("dir() should work on a registered class")
    };
    listing
        .iter()
        .filter_map(|entry| entry.extract::<String>().ok())
        .filter(|name| !name.starts_with('_') || CONTRACT_DUNDERS.contains(&name.as_str()))
        .collect()
}

fn describe(label: &str, extra: &BTreeSet<String>) -> String {
    format!("{label}: {extra:?}")
}

/// Records what `native` has that `declared` lacks, and the reverse.
fn compare(
    owner: &str,
    (extra, absent): (&str, &str),
    native: &BTreeSet<String>,
    declared: &BTreeSet<String>,
    problems: &mut Vec<String>,
) {
    let unlisted: BTreeSet<String> = native.difference(declared).cloned().collect();
    let stale: BTreeSet<String> = declared.difference(native).cloned().collect();
    if !unlisted.is_empty() {
        problems.push(describe(&format!("{owner}: {extra}"), &unlisted));
    }
    if !stale.is_empty() {
        problems.push(describe(&format!("{owner}: {absent}"), &stale));
    }
}

/// Compares every stubbed class's members with the registered class's.
fn compare_classes(
    owner: &str,
    space: &Bound<'_, PyAny>,
    stub: &stubs::Stub,
    problems: &mut Vec<String>,
) {
    for (class, expected) in &stub.members {
        let Ok(value) = space.getattr(class.as_str()) else {
            continue;
        };
        if value.is_instance_of::<PyType>() {
            compare(
                &format!("{owner}{class}"),
                ("members not in the stub", "stub members the class lacks"),
                &native_members(&value),
                expected,
                problems,
            );
        }
    }
}

#[cfg(feature = "full")]
#[test]
fn every_namespace_matches_its_stub_name_for_name_and_member_for_member() {
    let root = python_root();
    let native_stub = reexports(&read(&root.join("_native.pyi")));
    let root_stub = stubs::parse_root_stub(&root);
    let mut problems = Vec::new();
    let mut every_stub_name = BTreeSet::new();
    with_native(|module| {
        let spaces = namespaces(module);
        assert!(!spaces.is_empty(), "no namespaces were registered");
        for (name, space) in &spaces {
            let stub = parse_stub(&read(&root.join(name).join("__init__.pyi")));
            every_stub_name.extend(stub.names.iter().cloned());
            let declared: BTreeSet<String> = stub
                .names
                .difference(&stub.protocols)
                .filter(|each| !PYTHON_DEFINED.contains(&(name.as_str(), each.as_str())))
                .cloned()
                .collect();
            compare(
                name,
                (
                    "registered but not in its stub",
                    "in its stub but not registered",
                ),
                &non_module_names(space),
                &declared,
                &mut problems,
            );
            compare_classes(&format!("{name}."), space, &stub, &mut problems);
        }
        // A class registered on the root and re-exported by a namespace is
        // declared by that namespace's stub; the namespaces themselves are
        // registered as modules.
        let space_names: BTreeSet<String> = spaces.iter().map(|(name, _)| name.clone()).collect();
        let native = non_module_names(module.as_any());
        let named_by_the_root: BTreeSet<String> = native_stub
            .iter()
            .filter(|name| !space_names.contains(*name))
            .cloned()
            .collect();
        let unlisted: BTreeSet<String> = native
            .iter()
            .filter(|name| !named_by_the_root.contains(*name) && !every_stub_name.contains(*name))
            .cloned()
            .collect();
        let stale: BTreeSet<String> = named_by_the_root.difference(&native).cloned().collect();
        if !unlisted.is_empty() {
            problems.push(describe("root: registered but in no stub", &unlisted));
        }
        if !stale.is_empty() {
            problems.push(describe("root: in _native.pyi but not registered", &stale));
        }
        compare_classes("", module.as_any(), &root_stub, &mut problems);
    });
    assert!(
        problems.is_empty(),
        "stubs and extension disagree:\n{}",
        problems.join("\n")
    );
}

#[cfg(feature = "full")]
#[test]
fn every_namespace_lists_what_it_exports() {
    let root = python_root();
    let mut problems = Vec::new();
    with_native(|module| {
        for (name, space) in namespaces(module) {
            let exported = dunder_all(&read(&root.join(&name).join("__init__.py")));
            let missing: BTreeSet<String> = non_module_names(&space)
                .difference(&exported)
                .cloned()
                .collect();
            if !missing.is_empty() {
                problems.push(describe(
                    &format!("{name}: registered but not in __all__"),
                    &missing,
                ));
            }
        }
    });
    assert!(
        problems.is_empty(),
        "__all__ is incomplete:\n{}",
        problems.join("\n")
    );
}

#[cfg(feature = "full")]
#[test]
fn the_catalog_and_the_registered_operations_name_each_other() {
    let features: BTreeSet<String> =
        read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml"))
            .lines()
            .skip_while(|line| line.trim() != "[features]")
            .skip(1)
            .take_while(|line| !line.starts_with('['))
            .filter_map(|line| line.split_once('=').map(|(name, _)| name.trim().to_owned()))
            .collect();
    let mut problems = Vec::new();
    let mut seen = BTreeSet::new();
    let mut rows = BTreeSet::new();
    for capability in crate::catalog::capabilities() {
        if !seen.insert((capability.domain, capability.name)) {
            problems.push(format!(
                "{}.{} is catalogued twice",
                capability.domain, capability.name
            ));
        }
        if !features.contains(capability.feature) {
            problems.push(format!(
                "{}.{} names feature {:?}, which Cargo.toml does not declare",
                capability.domain, capability.name, capability.feature
            ));
        }
        rows.insert((capability.domain.to_owned(), capability.name.to_owned()));
    }
    let mut registered = BTreeSet::new();
    with_native(|module| {
        for (domain, space) in namespaces(module) {
            if UNCATALOGUED_NAMESPACES.contains(&domain.as_str()) {
                continue;
            }
            for (name, value) in public_attributes(&space) {
                if value.is_instance_of::<PyCFunction>() {
                    registered.insert((domain.clone(), name));
                }
            }
        }
        for capability in crate::catalog::capabilities().filter(|each| each.eager) {
            let present = module
                .getattr(capability.domain)
                .and_then(|space| space.getattr(capability.name))
                .is_ok();
            if !present {
                problems.push(format!(
                    "{}.{} is catalogued but absent",
                    capability.domain, capability.name
                ));
            }
        }
    });
    for missing in registered.difference(&rows) {
        problems.push(format!(
            "{}.{} is registered but has no catalog row",
            missing.0, missing.1
        ));
    }
    assert!(
        problems.is_empty(),
        "catalog and registration disagree:\n{}",
        problems.join("\n")
    );
}
