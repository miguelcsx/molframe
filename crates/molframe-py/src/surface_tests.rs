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
use std::fs;
use std::path::{Path, PathBuf};

#[path = "surface_tests/stubs.rs"]
mod stubs;
use stubs::{
    CONTRACT_DUNDERS, UNCATALOGUED_NAMESPACES, dunder_all, identifier, parse_stub, python_root,
    read, reexports, source_root,
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
    let root_stub = parse_stub(&read(&root.join("__init__.pyi")));
    let mut problems = Vec::new();
    let mut every_stub_name = BTreeSet::new();
    with_native(|module| {
        let spaces = namespaces(module);
        assert!(!spaces.is_empty(), "no namespaces were registered");
        for (name, space) in &spaces {
            let stub = parse_stub(&read(&root.join(name).join("__init__.pyi")));
            every_stub_name.extend(stub.names.iter().cloned());
            let declared: BTreeSet<String> =
                stub.names.difference(&stub.protocols).cloned().collect();
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

/// `#[pyfunction]`s that do not release the GIL, and why that is right: each
/// is a table lookup, a small expression constructor, or a single linear pass
/// over arrays already borrowed from `NumPy`, where handing the arrays across
/// the boundary would cost as much as the pass.
const HOLDS_THE_GIL: &[(&str, &str)] = &[
    ("bindings/geometry.rs", "centroid"),
    ("bindings/geometry.rs", "rmsd"),
    ("chemistry.rs", "element"),
    ("chemistry.rs", "vdw_radius"),
    ("chemistry/carbohydrates.rs", "snfg_symbol"),
    ("compare.rs", "tm_score"),
    ("compare.rs", "gdt_ts"),
    ("compare.rs", "gdt_ha"),
    ("crystal.rs", "assemblies"),
    ("crystal_reduction.rs", "reduce_cell"),
    ("query_aliases.rs", "complete"),
    ("selection_expr.rs", "(macro)"),
    ("selection_expr.rs", "chain"),
    ("selection_expr.rs", "residue"),
    ("selection_expr.rs", "atom"),
    ("selection_expr.rs", "within"),
    ("selection_expr.rs", "residues_within"),
];

fn rust_sources(directory: &Path, found: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(directory) else {
        panic!("{} should be readable", directory.display())
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            rust_sources(&path, found);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            found.push(path);
        }
    }
}

/// `(file, function)` for every `#[pyfunction]` whose body never releases the GIL.
fn functions_that_hold_the_gil() -> BTreeSet<(String, String)> {
    let source = source_root();
    let mut files = Vec::new();
    rust_sources(&source, &mut files);
    let mut holding = BTreeSet::new();
    for path in files {
        let relative = path
            .strip_prefix(&source)
            .map(|path| path.to_string_lossy().replace('\\', "/"))
            .unwrap_or_default();
        if relative.ends_with("_tests.rs") || relative.starts_with("native_source") {
            continue;
        }
        let text = read(&path);
        for (start, _) in text.match_indices("#[pyfunction]") {
            let tail = &text[start..];
            let Some(position) = tail.find("fn ") else {
                continue;
            };
            // A function written inside a macro has no name to list.
            let name = match identifier(&tail[position + 3..]) {
                name if name.is_empty() => "(macro)".to_owned(),
                name => name,
            };
            let end = tail[position..]
                .find("\n}\n")
                .map_or(tail.len(), |end| position + end);
            let body = &tail[..end];
            let releases = body.contains("detach(") || body.contains("run(");
            if !releases {
                holding.insert((relative.clone(), name));
            }
        }
    }
    holding
}

#[test]
fn every_pyfunction_releases_the_gil_unless_it_is_listed_as_trivial() {
    let holding = functions_that_hold_the_gil();
    let allowed: BTreeSet<(String, String)> = HOLDS_THE_GIL
        .iter()
        .map(|(file, name)| ((*file).to_owned(), (*name).to_owned()))
        .collect();
    let unlisted: Vec<_> = holding.difference(&allowed).collect();
    assert!(
        unlisted.is_empty(),
        "these #[pyfunction]s hold the GIL for the whole call; release it with \
         py.detach (or governed::run) or list them in HOLDS_THE_GIL with a reason: {unlisted:?}"
    );
    let stale: Vec<_> = allowed.difference(&holding).collect();
    assert!(
        stale.is_empty(),
        "HOLDS_THE_GIL lists functions that now release the GIL or no longer exist: {stale:?}"
    );
}

#[test]
fn no_python_source_exceeds_the_file_cap() {
    const CAP: usize = 500;
    let mut stack = vec![python_root()];
    let mut oversized = Vec::new();
    while let Some(directory) = stack.pop() {
        let Ok(entries) = fs::read_dir(&directory) else {
            panic!("{} should be readable", directory.display())
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().is_none_or(|name| name != "__pycache__") {
                    stack.push(path);
                }
            } else if path
                .extension()
                .is_some_and(|extension| extension == "py" || extension == "pyi")
            {
                let lines = read(&path).lines().count();
                if lines > CAP {
                    oversized.push(format!("{} has {lines} lines", path.display()));
                }
            }
        }
    }
    assert!(
        oversized.is_empty(),
        "over the {CAP}-line cap: {oversized:?}"
    );
}

#[test]
fn the_bindings_depend_on_the_facade_and_the_two_binding_libraries_only() {
    let manifest = read(&Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml"));
    let dependencies: BTreeSet<String> = manifest
        .lines()
        .skip_while(|line| line.trim() != "[dependencies]")
        .skip(1)
        .take_while(|line| !line.starts_with('['))
        .filter_map(|line| {
            line.split_once('=')
                .map(|(name, _)| name.trim().split('.').next().unwrap_or_default().to_owned())
        })
        .collect();
    let expected: BTreeSet<String> = ["molframe", "numpy", "pyo3"].map(str::to_owned).into();
    assert_eq!(
        dependencies, expected,
        "bindings stay thin: every kernel is reached through the facade crate"
    );
}
