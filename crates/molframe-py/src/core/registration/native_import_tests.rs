//! Every name a stub imports from the compiled extension is declared somewhere.
//!
//! The extension is one module in Rust and 400-odd stub files in Python, and a
//! stub reaches the extension by importing from `_native`. Nothing checked that
//! the name exists on the other side, so nine of them — `CoordinateEditor`,
//! `DictionaryFull`, `MissingResidue`, `FxEvaluation`, `MrcBrickError`,
//! `TrajectoryInterpolation`, `TrajectoryInterpolationError`,
//! `classify_ramachandran` and `AtomSelection` — were imported by a stub and
//! declared by none: a type checker reported the module as missing, and no test
//! noticed because no test read a `.pyi`.
//!
//! Only one direction is asserted, and it is the decidable one: everything a
//! stub imports must be declared. A name the extension binds and no stub names
//! is not a failure — a stub is free to leave out surface nobody types against.
//!
//! What this does not decide is *where* a name is imported from: a stub that
//! imports `AlgorithmId` from `.core` when only `.core.contract` declares it
//! passes, because the name exists. That residue is measured — 69 names in each
//! of the three facade stubs (`__init__.pyi`, `_facade_models.pyi`,
//! `_facade_structure_io.pyi`) and 39 singletons elsewhere, out of the 3 304 stub
//! import edges — and it is a stale copy of an earlier surface in a header those
//! three share, not a registration gap. Deciding it needs the module graph
//! resolved in Rust, the way `namespaces_tests` resolves a namespace's surface;
//! until then a stub-only edit here would be unverifiable, because the
//! repository has no type checker and no interpreter to run one in.

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

#[test]
fn every_name_a_stub_imports_from_the_extension_is_declared() {
    let root = package_root();
    let mut paths = Vec::new();
    stubs(&root, &mut paths);
    paths.sort();
    assert!(!paths.is_empty(), "no stub found under {}", root.display());

    let stubs: Vec<(PathBuf, String)> = paths
        .iter()
        .filter_map(|path| {
            std::fs::read_to_string(path)
                .ok()
                .map(|text| (path.clone(), text))
        })
        .collect();
    let declared: BTreeSet<String> = stubs
        .iter()
        .flat_map(|(_, text)| declarations(text))
        .collect();

    let mut failures = String::new();
    for (path, text) in &stubs {
        for (module, name) in extension_imports(text) {
            if !declared.contains(&name) {
                let _ = writeln!(
                    failures,
                    "{}: `{name}` is imported from `{module}` but no stub declares it",
                    path.display()
                );
            }
        }
    }
    assert!(
        failures.is_empty(),
        "the extension surface and the stubs disagree:\n{failures}"
    );
}

/// The `molframe` package, whose stubs describe the compiled extension.
fn package_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../python/molframe")
}

/// Every `.pyi` under `dir`, at any depth.
fn stubs(dir: &Path, found: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            stubs(&path, found);
        } else if path.extension().is_some_and(|extension| extension == "pyi") {
            found.push(path);
        }
    }
}

/// The `(module, name)` pairs a stub imports from the compiled extension.
fn extension_imports(text: &str) -> Vec<(String, String)> {
    let mut found = Vec::new();
    for statement in statements(text) {
        let Some(rest) = statement.strip_prefix("from ") else {
            continue;
        };
        let Some((module, list)) = rest.split_once(" import ") else {
            continue;
        };
        let module = module.trim();
        if !is_extension_module(module) {
            continue;
        }
        for name in imported_names(list) {
            found.push((module.to_owned(), name));
        }
    }
    found
}

/// Whether a `from …` module names the compiled extension, flat or as one of the
/// namespaces it exposes as a submodule.
fn is_extension_module(module: &str) -> bool {
    module
        .trim_start_matches('.')
        .split('.')
        .next()
        .is_some_and(|head| head == "_native")
}

/// The names an import list binds, with the parenthesised form unwrapped and any
/// `as` alias dropped. A star and the `__all__` re-export bind no name.
fn imported_names(list: &str) -> Vec<String> {
    list.trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split(',')
        .filter_map(|entry| {
            let name = entry.trim().split(" as ").next().unwrap_or_default().trim();
            (is_identifier(name) && name != "*" && name != "__all__").then(|| name.to_owned())
        })
        .collect()
}

/// The module's logical lines, so that an import list split over several lines
/// is read as the one statement it is.
fn statements(text: &str) -> Vec<String> {
    let mut statements = Vec::new();
    let mut open: Option<String> = None;
    for line in text.lines() {
        let joined = match open.take() {
            Some(mut buffer) => {
                buffer.push(' ');
                buffer.push_str(line.trim());
                buffer
            }
            None => line.trim().to_owned(),
        };
        if joined.contains("import (") && !joined.contains(')') {
            open = Some(joined);
        } else {
            statements.push(joined);
        }
    }
    statements.extend(open);
    statements
}

/// Every name the stub declares: a `class`, a `def`, an annotation, an
/// assignment, or an alias. A comment and an import declare nothing.
fn declarations(text: &str) -> BTreeSet<String> {
    text.lines().filter_map(declared_name).collect()
}

fn declared_name(line: &str) -> Option<String> {
    let trimmed = line.trim();
    if trimmed.starts_with('#') {
        return None;
    }
    for keyword in ["class ", "def "] {
        if let Some(rest) = trimmed.strip_prefix(keyword) {
            let name = rest.split(['(', ':', '[']).next().unwrap_or_default();
            return is_identifier(name).then(|| name.to_owned());
        }
    }
    let (name, rest) = trimmed.split_once([':', '='])?;
    let name = name.trim();
    (is_identifier(name) && !rest.is_empty()).then(|| name.to_owned())
}

fn is_identifier(name: &str) -> bool {
    !name.is_empty()
        && !name.starts_with(|character: char| character.is_ascii_digit())
        && name
            .chars()
            .all(|character| character.is_alphanumeric() || character == '_')
}
