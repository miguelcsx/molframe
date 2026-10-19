//! Gates on the Rust sources: the GIL discipline of every `#[pyfunction]`, the size of
//! the Python sources and the crate's dependency list.
//!
//! Cost: one read of every source file under `src/` and `python/`; no extension is
//! loaded.

use super::stubs::{identifier, python_root, read, source_root};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

/// `#[pyfunction]`s that do not release the GIL, and why that is right: each
/// is a table lookup, a handle constructor, a small expression constructor, or a
/// single linear pass over arrays already borrowed from `NumPy`, where handing
/// the arrays across the boundary would cost as much as the pass.
const HOLDS_THE_GIL: &[(&str, &str)] = &[
    ("bindings/geometry.rs", "centroid"),
    ("bindings/geometry.rs", "rmsd"),
    ("chemistry.rs", "element"),
    ("chemistry/roles.rs", "polymer_atom_roles"),
    ("chemistry.rs", "vdw_radius"),
    ("chemistry/carbohydrates.rs", "snfg_symbol"),
    ("crystal.rs", "assemblies"),
    ("crystal_assembly.rs", "assembly_chunks"),
    ("crystal_ncs.rs", "ncs_operators"),
    ("modelcif.rs", "pae"),
    ("modelcif.rs", "plddt_band_colors"),
    ("modelcif.rs", "plddt_band_names"),
    ("modelcif.rs", "plddt_bands"),
    ("interop/dlpack.rs", "coordinates"),
    ("crystal_reduction.rs", "reduce_cell"),
    ("query_aliases.rs", "complete"),
    ("selection_expr.rs", "(macro)"),
    ("selection_expr.rs", "chain"),
    ("selection_expr.rs", "residue"),
    ("selection_expr.rs", "element"),
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

fn is_identifier_character(character: char) -> bool {
    character.is_alphanumeric() || character == '_'
}

/// Every function in `text` with the source of its body, which ends at the
/// closing brace in the column the function's own line starts in.
fn function_bodies(text: &str) -> Vec<(String, &str)> {
    text.match_indices("fn ")
        .filter(|(start, _)| {
            text[..*start]
                .chars()
                .next_back()
                .is_none_or(|before| !is_identifier_character(before))
        })
        .map(|(start, _)| {
            let line_start = text[..start].rfind('\n').map_or(0, |newline| newline + 1);
            let indent = text[line_start..]
                .chars()
                .take_while(|character| *character == ' ')
                .count();
            let tail = &text[start..];
            let closing = format!("\n{}}}\n", " ".repeat(indent));
            let end = match tail.find(&closing) {
                Some(end) => end,
                None => tail.len(),
            };
            (identifier(&tail[3..]), &tail[..end])
        })
        .collect()
}

/// Whether `body` calls `name`.
fn calls(body: &str, name: &str) -> bool {
    let call = format!("{name}(");
    body.match_indices(&call).any(|(start, _)| {
        body[..start]
            .chars()
            .next_back()
            .is_none_or(|before| !is_identifier_character(before))
    })
}

/// The functions of one file that release the GIL, directly or by calling a
/// function of the same file that does.
fn releasing_functions(text: &str) -> BTreeSet<String> {
    let functions = function_bodies(text);
    let mut releasing: BTreeSet<String> = functions
        .iter()
        .filter(|(_, body)| body.contains("detach(") || body.contains("run("))
        .map(|(name, _)| name.clone())
        .collect();
    loop {
        let before = releasing.len();
        for (name, body) in &functions {
            if !releasing.contains(name) && releasing.iter().any(|callee| calls(body, callee)) {
                releasing.insert(name.clone());
            }
        }
        if releasing.len() == before {
            return releasing;
        }
    }
}

/// `(file, function)` for every `#[pyfunction]` that never releases the GIL,
/// directly or through a helper in its own file.
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
        let releasing = releasing_functions(&text);
        for (start, _) in text.match_indices("#[pyfunction]") {
            let tail = &text[start..];
            let Some(position) = tail.find("fn ") else {
                continue;
            };
            let name = identifier(&tail[position + 3..]);
            let end = tail[position..]
                .find("\n}\n")
                .map_or(tail.len(), |end| position + end);
            let body = &tail[..end];
            let releases = if name.is_empty() {
                // A function written inside a macro has no name to look up.
                body.contains("detach(") || body.contains("run(")
            } else {
                releasing.contains(&name)
            };
            if !releases {
                let listed = if name.is_empty() {
                    "(macro)".to_owned()
                } else {
                    name
                };
                holding.insert((relative.clone(), listed));
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
