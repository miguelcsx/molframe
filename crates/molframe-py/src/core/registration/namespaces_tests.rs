//! The registered namespaces, checked against the Python surface they back.
//!
//! Each namespace registers a list of names and each `python/molframe/<ns>/`
//! module binds them, and the two lists are written by different hands in
//! different languages with nothing in between. Before this check, `core` was
//! registered with 153 names while its stub declared 160, and five of the
//! metadata types the module re-exported reached no stub at all.
//!
//! Only the direction that can be decided is asserted: a registered name must
//! be reachable. The converse does not hold — `molframe.core` deliberately
//! exposes more than the `core` namespace registers — so it is not checked.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use super::NAMESPACES;
use super::namespace_overlays;

#[test]
fn every_registered_name_is_declared_by_the_namespace_that_carries_it() {
    let root = package_root();
    let mut failures = String::new();

    for (namespace, exports) in NAMESPACES {
        // What the extension publishes is not `exports` alone. `register` adds
        // the namespace's overlays, aliases and child modules and lists all four
        // in `__all__`, so a name from any of them is reachable at the namespace
        // path and the stub has to name it too. Checking `exports` alone is what
        // left every overlay name — `BasePair`, `CationPi`, `Contact` — declared
        // on no stub at all while the runtime bound them.
        let surface: Vec<&str> = exports
            .iter()
            .chain(namespace_overlays::for_namespace(namespace))
            .chain(namespace_overlays::alias_names(namespace))
            .chain(namespace_overlays::children(namespace))
            .copied()
            .collect();
        // A module that re-exports the native namespace wholesale has no list of
        // its own, and the runtime has the names because the extension binds
        // them. The stub has to declare them somewhere, and there are two places
        // it may: the namespace's own stub, which gives the names real
        // signatures, and the stub of the extension submodule it defers to. A
        // namespace that defers must be declared by one or the other.
        let wholesale = format!("from .._native.{namespace} import *");
        for extension in ["py", "pyi"] {
            let file = module_path(&root, namespace, extension);
            let text = match std::fs::read_to_string(&file) {
                Ok(text) => text,
                Err(error) => {
                    let _ = writeln!(failures, "{}: {error}", file.display());
                    continue;
                }
            };
            let defers = text.contains(&wholesale);
            if defers && extension == "py" {
                // The runtime module gets the names from the extension itself.
                continue;
            }
            let mut declared = bindings(&text).to_owned();
            if defers {
                let carrier = carrier_stub(&root, namespace);
                match std::fs::read_to_string(&carrier) {
                    Ok(text) => {
                        declared.push('\n');
                        declared.push_str(bindings(&text));
                    }
                    Err(error) => {
                        let _ = writeln!(failures, "{}: {error}", carrier.display());
                        continue;
                    }
                }
            }
            for export in surface.iter().copied() {
                if !mentions(&declared, export) {
                    let _ = writeln!(
                        failures,
                        "{namespace}: `{export}` is registered but never named in {}",
                        file.display()
                    );
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "the registration and the Python surface disagree:\n{failures}"
    );
}

/// The `molframe` package, whose namespace modules this file's list backs.
fn package_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../python/molframe")
}

fn module_path(root: &Path, namespace: &str, extension: &str) -> PathBuf {
    root.join(namespace).join(format!("__init__.{extension}"))
}

/// The stub of the extension submodule that carries a namespace's names.
fn carrier_stub(root: &Path, namespace: &str) -> PathBuf {
    root.join("_native").join(format!("{namespace}.pyi"))
}

/// The part of a module that establishes its bindings, which is every part
/// before it restates them. Only a declaration or an import binds, so a name
/// that survives solely as a string in `__all__` is a name the module no longer
/// has — and would raise on import.
///
/// The cut is at the assignment, not at the first mention. A module that
/// re-exports the native `__all__` names it near the top, and cutting there
/// would hide every declaration below it: the whole of `validate/__init__.pyi`
/// was invisible to this check for exactly that reason.
fn bindings(text: &str) -> &str {
    let mut cut = text.len();
    let mut start = 0;
    for line in text.split('\n') {
        let restated = line.trim_start();
        if restated.starts_with("__all__") && restated.contains(['=', ':']) {
            cut = start;
            break;
        }
        start += line.len() + 1;
    }
    &text[..cut]
}

/// Whether `name` occurs in `text` as a whole identifier, so that `Column` is
/// not satisfied by `ColumnKind` and `read` is not satisfied by `read_mmtf`.
fn mentions(text: &str, name: &str) -> bool {
    let mut from = 0;
    while let Some(offset) = text[from..].find(name) {
        let start = from + offset;
        let end = start + name.len();
        if !is_identifier_char(text[..start].chars().next_back())
            && !is_identifier_char(text[end..].chars().next())
        {
            return true;
        }
        from = end;
    }
    false
}

fn is_identifier_char(character: Option<char>) -> bool {
    match character {
        Some(character) => character.is_alphanumeric() || character == '_',
        None => false,
    }
}
