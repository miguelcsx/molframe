//! SASA as a method, and the ensemble gather it shares a guard with.

use super::*;
use crate::prelude::*;

const DENSE: &str = "\
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1       5.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       6.458   0.000   0.000  1.00 10.00           C
ENDMDL
END
";

const RAGGED: &str = "\
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1       5.000   0.000   0.000  1.00 10.00           N
ENDMDL
END
";

fn read(text: &str) -> crate::Structure {
    match crate::read_bytes(text.as_bytes().to_vec(), Some("e.pdb"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture read failed: {findings:?}"),
    }
}

#[test]
fn sasa_gives_one_positive_area_per_atom_and_sums_per_residue() {
    let structure = read(DENSE);
    let areas = structure.sasa().expect("sasa");
    assert_eq!(areas.len(), 2);
    assert!(areas.iter().all(|area| *area > 0.0));
    let residues = structure.residue_sasa().expect("residue sasa");
    assert_eq!(residues.len(), 1);
    assert!((residues[0] - areas.iter().sum::<f64>()).abs() < 1e-9);
}

#[test]
fn sasa_measures_the_requested_model() {
    let structure = read(DENSE);
    let options = SasaOptions {
        model: ModelIndex::new(1),
        ..SasaOptions::default()
    };
    let second = structure
        .sasa_with(&options, &ExecutionContext::default())
        .expect("second model");
    let first = structure.sasa().expect("first model");
    // The models are translations of each other, so the areas agree.
    assert!((first[0] - second[0]).abs() < 1e-6);
    let missing = SasaOptions {
        model: ModelIndex::new(9),
        ..SasaOptions::default()
    };
    assert!(matches!(
        structure.sasa_with(&missing, &ExecutionContext::default()),
        Err(SurfaceExtError::Structure(_))
    ));
}

#[test]
fn a_ragged_ensemble_is_refused_by_name() {
    let structure = read(RAGGED);
    let Err(SurfaceExtError::Structure(finding)) = structure.sasa() else {
        panic!("ragged ensemble must be refused");
    };
    assert!(finding.to_string().contains("ragged_models"));
    let Err(finding) = structure.ensemble_coordinates(None) else {
        panic!("ragged ensemble has no shared atom axis");
    };
    assert!(finding.to_string().contains("ragged_models"));
}

#[test]
fn a_dense_ensemble_gathers_into_one_block() {
    let structure = read(DENSE);
    let block = structure.ensemble_coordinates(None).expect("block");
    assert_eq!((block.models(), block.atoms()), (2, 2));
    assert_eq!(block.as_flat().len(), 12);
    assert_eq!(
        block.model(1).expect("model 1")[0].map(f32::to_bits),
        [5.0_f32, 0.0, 0.0].map(f32::to_bits)
    );
    let one = structure
        .ensemble_coordinates(Some(&[ModelIndex::new(1)]))
        .expect("subset");
    assert_eq!((one.models(), one.atoms()), (1, 2));
}

#[test]
fn atoms_can_be_iterated_directly() {
    let structure = read(DENSE);
    let names: Vec<_> = structure
        .atoms()
        .into_iter()
        .filter_map(AtomRef::name)
        .collect();
    assert_eq!(names, ["N", "CA"]);
}
