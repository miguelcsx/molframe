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

/// The protein-only alpha (A) and beta (B) chains of 4HHB; chain B is moved by
/// `shift` Angstrom along x when it is not zero.
pub fn haemoglobin_dimer(directory: &Path, name: &str, shift: f64) -> PathBuf {
    let source = std::fs::read_to_string(bench_data("4hhb.pdb")).expect("bundled 4HHB");
    let mut text = String::new();
    for line in source.lines() {
        if !line.starts_with("ATOM") || !matches!(line.as_bytes()[21], b'A' | b'B') {
            continue;
        }
        let mut x: f64 = line[30..38].trim().parse().expect("x column");
        if line.as_bytes()[21] == b'B' {
            x += shift;
        }
        let _ = writeln!(text, "{}{x:8.3}{}", &line[..30], &line[38..]);
    }
    text.push_str("END\n");
    let path = directory.join(name);
    std::fs::write(&path, text).expect("write dimer");
    path
}

/// The peptide rotated a quarter turn about z and moved, as a rigid body.
pub fn rigidly_moved_peptide(directory: &Path) -> PathBuf {
    let source = std::fs::read_to_string(hydrogenated_peptide(directory)).expect("peptide");
    let mut text = String::new();
    for line in source.lines() {
        if line.starts_with("ATOM") {
            let x: f64 = line[30..38].trim().parse().expect("x");
            let y: f64 = line[38..46].trim().parse().expect("y");
            let _ = writeln!(
                text,
                "{}{:8.3}{:8.3}{}",
                &line[..30],
                -y + 5.0,
                x - 3.0,
                &line[46..]
            );
        } else {
            text.push_str(line);
            text.push('\n');
        }
    }
    let path = directory.join("moved.pdb");
    std::fs::write(&path, text).expect("write moved peptide");
    path
}

/// The peptide with the two ring-flip-equivalent halves of its first tyrosine
/// swapped, atom names left in place, so only a symmetry-aware RMSD sees through it.
pub fn ring_flipped_peptide(directory: &Path) -> (PathBuf, String) {
    let source = std::fs::read_to_string(hydrogenated_peptide(directory)).expect("peptide");
    let lines: Vec<&str> = source.lines().collect();
    let first = lines
        .iter()
        .position(|line| line.starts_with("ATOM") && &line[17..20] == "TYR")
        .expect("a tyrosine");
    let sequence = lines[first][22..26].to_owned();
    let residue: Vec<usize> = (first..lines.len())
        .take_while(|&index| lines[index].starts_with("ATOM") && lines[index][22..26] == sequence)
        .collect();
    let find = |name: &str| {
        *residue
            .iter()
            .find(|&&index| lines[index][12..16].trim() == name)
            .unwrap_or_else(|| panic!("tyrosine atom {name}"))
    };
    let mut out: Vec<String> = lines.iter().map(|line| (*line).to_owned()).collect();
    for (left, right) in [
        ("CD1", "CD2"),
        ("CE1", "CE2"),
        ("HD1", "HD2"),
        ("HE1", "HE2"),
    ] {
        let (a, b) = (find(left), find(right));
        out[a] = format!(
            "{}{}{}",
            &lines[a][..30],
            &lines[b][30..54],
            &lines[a][54..]
        );
        out[b] = format!(
            "{}{}{}",
            &lines[b][..30],
            &lines[a][30..54],
            &lines[b][54..]
        );
    }
    let path = directory.join("flipped.pdb");
    std::fs::write(&path, out.join("\n") + "\n").expect("write flipped peptide");
    (path, sequence.trim().to_owned())
}
