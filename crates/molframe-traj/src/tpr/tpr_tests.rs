use super::*;
use std::path::Path;

#[test]
fn rejects_non_tpr_input() {
    assert!(matches!(
        parse_tpr(b"not a tpr"),
        Err(TprError::Truncated { .. })
    ));
}

#[test]
fn reads_modern_single_atom_topology_without_defaults() {
    let topology = match parse_tpr(&single_atom_tpr()) {
        Ok(topology) => topology,
        Err(error) => panic!("synthetic TPR failed: {error}"),
    };
    assert_eq!(topology.header.format_version, 138);
    assert_eq!(topology.header.generation, 29);
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(topology.residues.len(), 1);
    assert_eq!(topology.atoms[0].name, "CA");
    assert_eq!(topology.atoms[0].atom_type, "CT");
    assert_eq!(topology.atoms[0].atomic_number, Some(6));
    assert!((topology.atoms[0].mass - 12.011).abs() < 1.0e-5);
    assert!((topology.atoms[0].charge + 0.1).abs() < 1.0e-5);
    assert_eq!(topology.residues[0].name, "GLY");
    assert_eq!(topology.residues[0].molecule_type, "Protein_A");
    assert!(topology.bonds.is_empty());
}

#[test]
fn reads_external_version_corpus_when_configured() {
    let Ok(root) = std::env::var("MOLFRAME_TPR_CORPUS") else {
        return;
    };
    let entries = std::fs::read_dir(&root).expect("TPR corpus directory");
    let mut parsed = 0;
    for entry in entries {
        let path = entry.expect("corpus entry").path();
        if path.is_dir() {
            parsed += parse_directory(&path);
        } else if path.extension().is_some_and(|extension| extension == "tpr") {
            parsed += usize::from(verify_file(&path));
        }
    }
    assert!(parsed > 0, "TPR corpus contained no files");
}

fn parse_directory(root: &Path) -> usize {
    let mut parsed = 0;
    for entry in std::fs::read_dir(root).expect("nested TPR corpus directory") {
        let path = entry.expect("nested corpus entry").path();
        if path.is_dir() {
            parsed += parse_directory(&path);
        } else if path.extension().is_some_and(|extension| extension == "tpr") {
            parsed += usize::from(verify_file(&path));
        }
    }
    parsed
}

fn verify_file(path: &Path) -> bool {
    let bytes = std::fs::read(path).expect("TPR fixture bytes");
    let topology = match parse_tpr(&bytes) {
        Ok(topology) => topology,
        Err(TprError::UnsupportedVersion(_) | TprError::UnsupportedBetaSerializer) => return false,
        Err(error) => panic!("{} failed: {error}", path.display()),
    };
    assert_eq!(topology.atoms.len(), topology.header.atom_count);
    assert!(!topology.residues.is_empty());
    assert!(
        topology.bonds.iter().all(|bond| {
            bond.atom_a < topology.atoms.len() && bond.atom_b < topology.atoms.len()
        })
    );
    true
}

fn single_atom_tpr() -> Vec<u8> {
    let mut bytes = Vec::new();
    classic_string(&mut bytes, "VERSION synthetic");
    int(&mut bytes, 4);
    int(&mut bytes, 138);
    int(&mut bytes, 29);
    classic_string(&mut bytes, "release");
    for value in [1, 0] {
        int(&mut bytes, value);
    }
    int(&mut bytes, 0);
    real(&mut bytes, 0.0);
    for value in [0, 1, 0, 0, 0, 0] {
        int(&mut bytes, value);
    }
    bytes.extend_from_slice(&0_i64.to_be_bytes());

    let symbols = ["system", "Protein_A", "CA", "CT", "GLY"];
    int(&mut bytes, 5);
    for symbol in symbols {
        modern_string(&mut bytes, symbol);
    }
    int(&mut bytes, 0);
    int(&mut bytes, 1);
    int(&mut bytes, 0);
    bytes.extend_from_slice(&12.0_f64.to_be_bytes());
    real(&mut bytes, 1.0);
    int(&mut bytes, 1);
    int(&mut bytes, 1);
    int(&mut bytes, 1);
    int(&mut bytes, 1);
    real(&mut bytes, 12.011);
    real(&mut bytes, -0.1);
    real(&mut bytes, 12.011);
    real(&mut bytes, -0.1);
    bytes.extend_from_slice(&0_u16.to_be_bytes());
    bytes.extend_from_slice(&0_u16.to_be_bytes());
    int(&mut bytes, 0);
    int(&mut bytes, 0);
    int(&mut bytes, 6);
    for index in [2, 3, 3, 4] {
        int(&mut bytes, index);
    }
    int(&mut bytes, 1);
    bytes.push(0);
    for _ in 0..95 {
        int(&mut bytes, 0);
    }
    int(&mut bytes, 0);
    int(&mut bytes, 0);
    int(&mut bytes, 0);
    int(&mut bytes, 0);
    int(&mut bytes, 0);
    int(&mut bytes, 1);
    for value in [0, 1, 1, 0, 0] {
        int(&mut bytes, value);
    }
    bytes
}

fn int(bytes: &mut Vec<u8>, value: i32) {
    bytes.extend_from_slice(&value.to_be_bytes());
}

fn real(bytes: &mut Vec<u8>, value: f32) {
    bytes.extend_from_slice(&value.to_be_bytes());
}

fn classic_string(bytes: &mut Vec<u8>, value: &str) {
    let Ok(length) = i32::try_from(value.len()) else {
        panic!("test string length exceeds i32");
    };
    int(bytes, length);
    int(bytes, length);
    bytes.extend_from_slice(value.as_bytes());
    bytes.resize(bytes.len().next_multiple_of(4), 0);
}

fn modern_string(bytes: &mut Vec<u8>, value: &str) {
    let Ok(length) = u64::try_from(value.len()) else {
        panic!("test string length exceeds u64");
    };
    bytes.extend_from_slice(&length.to_be_bytes());
    bytes.extend_from_slice(value.as_bytes());
}
