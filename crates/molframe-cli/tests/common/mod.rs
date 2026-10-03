//! Shared harness for end-to-end command-line tests.
//!
//! Each test runs the built binary in a private directory against bundled real
//! data and reads the rows it prints, so what is asserted is what a user sees.

#![allow(dead_code)]

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn bench_data(name: &str) -> PathBuf {
    root().join("../molframe-bench/data").join(name)
}

pub fn amino_acid_ccd() -> PathBuf {
    root().join("../molframe-chem/data/CCD-amino-acids.cif")
}

/// Runs the binary with `arguments` in `directory` and returns its output.
pub fn run(directory: &Path, arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_molframe"))
        .current_dir(directory)
        .env("XDG_CONFIG_HOME", directory)
        .args(arguments)
        .output()
        .expect("run the CLI")
}

/// Runs and requires success, returning the tab-separated rows (header first).
pub fn rows(directory: &Path, arguments: &[&str]) -> Vec<Vec<String>> {
    let mut all = vec!["--format", "tsv"];
    all.splice(0..0, arguments.iter().copied());
    let output = run(directory, &all);
    assert!(
        output.status.success(),
        "{arguments:?} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .expect("utf-8 output")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

/// Model 1 of the bundled 2M7C NMR structure with its explicit hydrogens,
/// without the one non-standard residue that the twenty-residue dictionary
/// does not define.
pub fn hydrogenated_peptide(directory: &Path) -> PathBuf {
    let source = std::fs::read_to_string(bench_data("2m7c.pdb")).expect("bundled 2M7C");
    let mut text = String::new();
    for line in source.lines() {
        if line.starts_with("ENDMDL") {
            break;
        }
        if line.starts_with("ATOM") {
            text.push_str(line);
            text.push('\n');
        }
    }
    text.push_str("END\n");
    let path = directory.join("peptide.pdb");
    std::fs::write(&path, text).expect("write peptide");
    path
}

/// The columns of `table`'s header, which must be the first row.
pub fn header(table: &[Vec<String>]) -> &[String] {
    &table[0]
}

pub fn column(table: &[Vec<String>], name: &str) -> usize {
    table[0]
        .iter()
        .position(|candidate| candidate == name)
        .unwrap_or_else(|| panic!("no column {name} in {:?}", table[0]))
}

/// The same peptide reflected through the yz plane: every L stereocentre
/// becomes a D one while all distances and angles are unchanged.
pub fn mirrored_peptide(directory: &Path) -> PathBuf {
    let source = std::fs::read_to_string(hydrogenated_peptide(directory)).expect("peptide");
    let mut text = String::new();
    for line in source.lines() {
        if line.starts_with("ATOM") {
            let x: f64 = line[30..38].trim().parse().expect("x column");
            let _ = writeln!(text, "{}{:8.3}{}", &line[..30], -x, &line[38..]);
        } else {
            text.push_str(line);
            text.push('\n');
        }
    }
    let path = directory.join("mirrored.pdb");
    std::fs::write(&path, text).expect("write mirrored peptide");
    path
}

/// The conventional N/CA/C/O heavy-backbone roles for amino acids.
pub fn backbone_roles(directory: &Path) -> PathBuf {
    let path = directory.join("roles.toml");
    std::fs::write(
        &path,
        "profile_id = \"amino-acid-backbone-test\"\n\
         [[rules]]\natom_name = \"N\"\ncomponent_kind = 1\nrole = 1\n\
         [[rules]]\natom_name = \"CA\"\ncomponent_kind = 1\nrole = 2\n\
         [[rules]]\natom_name = \"C\"\ncomponent_kind = 1\nrole = 4\n\
         [[rules]]\natom_name = \"O\"\ncomponent_kind = 1\nrole = 8\n",
    )
    .expect("write role profile");
    path
}

/// Number of ATOM records in a PDB file.
pub fn atom_records(path: &Path) -> usize {
    std::fs::read_to_string(path)
        .expect("read structure")
        .lines()
        .filter(|line| line.starts_with("ATOM"))
        .count()
}

/// Runs and returns the exit code with the rows printed, without requiring success.
pub fn rows_with_status(directory: &Path, arguments: &[&str]) -> (i32, Vec<Vec<String>>) {
    let mut all: Vec<&str> = arguments.to_vec();
    all.extend(["--format", "tsv"]);
    let output = run(directory, &all);
    let rows = String::from_utf8(output.stdout)
        .expect("utf-8 output")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect();
    (output.status.code().expect("exit code"), rows)
}
