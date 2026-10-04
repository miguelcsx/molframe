//! Reading the hand-written Python stubs and package sources.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

/// Special methods that are part of a class's contract and so appear in stubs.
pub(super) const CONTRACT_DUNDERS: &[&str] = &[
    "__and__",
    "__arrow_c_stream__",
    "__bool__",
    "__contains__",
    "__dlpack__",
    "__dlpack_device__",
    "__enter__",
    "__exit__",
    "__float__",
    "__getitem__",
    "__int__",
    "__invert__",
    "__iter__",
    "__len__",
    "__next__",
    "__or__",
    "__sub__",
];

/// Stub members that say nothing about the class's callable surface.
pub(super) const STUB_ONLY_MEMBERS: &[&str] = &[
    "__init__", "__new__", "__repr__", "__str__", "__hash__", "__eq__",
];

/// Namespaces whose functions are expression constructors, not operations, and
/// so are not catalogued: each builds a small selection expression.
pub(super) const UNCATALOGUED_NAMESPACES: &[&str] = &["sel", "query"];

pub(super) fn python_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/molframe")
}

pub(super) fn source_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("src")
}

pub(super) fn read(path: &Path) -> String {
    match fs::read_to_string(path) {
        Ok(text) => text,
        Err(error) => panic!("{} should be readable: {error}", path.display()),
    }
}

pub(super) fn identifier(text: &str) -> String {
    text.chars()
        .take_while(|character| character.is_alphanumeric() || *character == '_')
        .collect()
}

#[derive(Default)]
pub(super) struct Stub {
    pub(super) names: BTreeSet<String>,
    pub(super) protocols: BTreeSet<String>,
    pub(super) members: BTreeMap<String, BTreeSet<String>>,
}

/// The name a class-body line declares: a `def`, or an annotated attribute.
fn member_name(body: &str) -> Option<String> {
    if let Some(rest) = body.strip_prefix("def ") {
        return Some(identifier(rest));
    }
    let name = identifier(body);
    body[name.len()..]
        .trim_start()
        .starts_with(':')
        .then_some(name)
}

/// The public names and class members a stub file declares.
pub(super) fn parse_stub(text: &str) -> Stub {
    let mut stub = Stub::default();
    let mut current: Option<String> = None;
    let mut in_docstring = false;
    for line in text.lines() {
        // Prose inside a docstring declares nothing, whatever it looks like.
        let delimiters = line.matches("\"\"\"").count();
        if in_docstring || delimiters % 2 == 1 {
            if delimiters % 2 == 1 {
                in_docstring = !in_docstring;
            }
            continue;
        }
        if let Some(rest) = line.strip_prefix("class ") {
            let name = identifier(rest);
            if name.starts_with('_') {
                current = None;
                continue;
            }
            if rest.contains("Protocol") {
                stub.protocols.insert(name.clone());
            }
            stub.names.insert(name.clone());
            stub.members.entry(name.clone()).or_default();
            current = Some(name);
            continue;
        }
        if !line.starts_with(' ') && !line.starts_with(')') && !line.is_empty() {
            current = None;
            let declared = line
                .strip_prefix("def ")
                .map(identifier)
                .or_else(|| {
                    let name = identifier(line);
                    let after = line[name.len()..].trim_start();
                    (after.starts_with(':') || after.starts_with('=')).then_some(name)
                })
                .filter(|name| !name.is_empty() && !name.starts_with('_'));
            if let Some(name) = declared {
                stub.names.insert(name);
            }
            continue;
        }
        let (Some(class), Some(body)) = (&current, line.strip_prefix("    ")) else {
            continue;
        };
        if body.starts_with(' ') {
            continue;
        }
        let Some(member) = member_name(body) else {
            continue;
        };
        let contract = CONTRACT_DUNDERS.contains(&member.as_str());
        if !member.is_empty()
            && (!member.starts_with('_') || contract)
            && !STUB_ONLY_MEMBERS.contains(&member.as_str())
        {
            stub.members
                .entry(class.clone())
                .or_default()
                .insert(member);
        }
    }
    stub
}

/// The names `from . import (X as X, ...)` re-exports.
pub(super) fn reexports(text: &str) -> BTreeSet<String> {
    text.lines()
        .filter_map(|line| {
            let (left, right) = line.trim().trim_end_matches(',').split_once(" as ")?;
            (left == right && !left.starts_with('_')).then(|| left.to_owned())
        })
        .collect()
}

/// The names a package's `__all__` lists.
pub(super) fn dunder_all(text: &str) -> BTreeSet<String> {
    let Some(start) = text.find("__all__") else {
        return BTreeSet::new();
    };
    let tail = &text[start..];
    let Some(open) = tail.find('[') else {
        return BTreeSet::new();
    };
    let Some(close) = tail.find(']') else {
        return BTreeSet::new();
    };
    tail[open..close]
        .split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_owned)
        .collect()
}

/// The root package's stub: `__init__.pyi` and the private parts it imports from
/// (`_structure.pyi`, ...), which exist only to keep each file under the cap.
pub(super) fn parse_root_stub(root: &Path) -> Stub {
    let mut merged = parse_stub(&read(&root.join("__init__.pyi")));
    let Ok(entries) = fs::read_dir(root) else {
        panic!("{} should be readable", root.display())
    };
    let mut parts: Vec<PathBuf> = entries
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension().is_some_and(|extension| extension == "pyi")
                && path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with('_') && name != "_native.pyi")
        })
        .collect();
    parts.sort();
    for part in parts {
        let stub = parse_stub(&read(&part));
        merged.names.extend(stub.names);
        merged.protocols.extend(stub.protocols);
        merged.members.extend(stub.members);
    }
    merged
}
