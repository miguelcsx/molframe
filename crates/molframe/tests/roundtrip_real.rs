//! Round-trip fidelity of real deposited files across cif, bcif and pdb:
//! text and binary mmCIF must be lossless, PDB must agree within its column
//! precision, and every PDB loss must be reported rather than silent.

use molframe::formats::pdb::{PdbIdentifierNamespace, PdbOptions};
use molframe::{ReadOptions, Structure, StructureDifferenceOptions};
use std::path::PathBuf;

const PDB_TOLERANCE: f32 = 0.0005;

fn bench(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../molframe-bench/data")
        .join(name)
}

fn read(bytes: Vec<u8>, name: &str) -> Structure {
    match molframe::read_bytes(bytes, Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(failure) => panic!("{name}: {failure:?}"),
    }
}

fn load(name: &str) -> Structure {
    let bytes = std::fs::read(bench(name)).unwrap();
    read(bytes, name)
}

fn auth_options() -> PdbOptions {
    PdbOptions::new().namespace(PdbIdentifierNamespace::Auth)
}

fn difference(
    left: &Structure,
    right: &Structure,
    tolerance: f32,
) -> molframe::StructureDifference {
    molframe::structure_difference(
        left,
        right,
        StructureDifferenceOptions {
            coordinate_tolerance: tolerance,
        },
    )
    .unwrap()
}

fn identical(left: &Structure, right: &Structure, tolerance: f32) -> bool {
    difference(left, right, tolerance).is_empty()
}

#[test]
fn text_mmcif_rewrites_without_any_difference() {
    for name in ["1crn.cif", "1ubq.cif", "4hhb.cif"] {
        let source = load(name);
        let text = molframe::write_mmcif(&source).unwrap();
        let again = read(text.into_bytes(), "again.cif");
        assert!(identical(&source, &again, 0.0), "{name}");
    }
}

#[test]
fn binary_mmcif_matches_text_mmcif_exactly() {
    for name in ["1crn.cif", "1ubq.cif", "4hhb.cif", "1aon.bcif"] {
        let source = load(name);
        let text = molframe::write_mmcif(&source).unwrap();
        let binary = read(molframe::write_bcif(&source).unwrap(), "again.bcif");
        assert!(identical(&source, &binary, 0.0), "{name}");
        assert_eq!(molframe::write_mmcif(&binary).unwrap(), text, "{name}");
    }
}

#[test]
fn pdb_projection_of_mmcif_agrees_within_column_precision() {
    for name in ["1crn.cif", "1ubq.cif", "4hhb.cif"] {
        let source = load(name);
        let pdb = molframe::write_pdb(&source, &auth_options()).unwrap();
        let again = read(pdb.into_bytes(), "again.pdb");
        let diff = difference(&source, &again, PDB_TOLERANCE);
        assert_eq!(atom_fields(&source), atom_fields(&again), "{name}");
        assert_eq!(diff.changed_positions, 0, "{name}");
    }
}

#[test]
fn pdb_rewrites_to_identical_coordinate_records() {
    for name in ["1crn.pdb", "1ubq.pdb", "4hhb.pdb", "2m7c.pdb"] {
        let source = load(name);
        let first = molframe::write_pdb(&source, &auth_options()).unwrap();
        let again = read(first.clone().into_bytes(), "again.pdb");
        let second = molframe::write_pdb(&again, &auth_options()).unwrap();
        assert_eq!(
            coordinate_lines(&first),
            coordinate_lines(&second),
            "{name}"
        );
    }
}

#[test]
fn pdb_to_mmcif_keeps_the_chain_count() {
    for name in ["1ubq.pdb", "4hhb.pdb", "1crn.pdb"] {
        let source = load(name);
        let text = molframe::write_mmcif(&source).unwrap();
        let again = read(text.into_bytes(), "again.cif");
        assert_eq!(
            source.engine().data().chains().count(),
            again.engine().data().chains().count(),
            "{name}"
        );
    }
}

#[test]
fn entity_bearing_mmcif_reports_its_pdb_loss() {
    let source = load("4hhb.cif");
    let (_, findings) = molframe::write_pdb_with_findings(&source, &auth_options()).unwrap();
    assert!(
        findings
            .iter()
            .any(|f| f.code().to_string().contains("W4101"))
    );
}

fn coordinate_lines(pdb: &str) -> Vec<&str> {
    pdb.lines()
        .filter(|l| {
            ["ATOM", "HETAT", "TER", "ANISO"]
                .iter()
                .any(|p| l.starts_with(p))
        })
        .collect()
}

/// Per-atom fields a PDB atom record carries; label chain ids are excluded
/// because PDB writes one namespace and reports the other as lost.
fn atom_fields(structure: &Structure) -> Vec<String> {
    (0..structure.atom_count())
        .map(|i| {
            let atom = structure
                .engine()
                .atom(molframe::AtomIndex::new(i))
                .unwrap();
            format!(
                "{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
                atom.auth_name(),
                atom.component_name(),
                atom.element(),
                atom.alt_label(),
                atom.b_factor(),
                atom.occupancy(),
                atom.formal_charge()
            )
        })
        .collect()
}
