//! Perception on multi-model reads, dense and ragged.
#![cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]

use molframe::{ReadOptions, Structure, read_bytes};

fn model(number: u32, b_factor: &str) -> String {
    format!(
        "MODEL     {number:>4}\n\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00{b_factor}           N\n\
ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00{b_factor}           C\n\
TER\n\
ATOM      3  N   GLY B   1      20.000   0.000   0.000  1.00{b_factor}           N\n\
ATOM      4  CA  GLY B   1      21.450   0.000   0.000  1.00{b_factor}           C\n\
TER\n\
ENDMDL\n"
    )
}

fn read(text: String) -> Structure {
    let Ok((structure, _)) = read_bytes(text.into_bytes(), Some("m.pdb"), &ReadOptions::new())
    else {
        panic!("multi-model fixture must read")
    };
    structure
}

#[test]
fn identical_models_with_two_chains_read_densely() {
    let structure = read(format!("{}{}END\n", model(1, " 10.00"), model(2, " 10.00")));
    assert_eq!(structure.model_count(), 2);
    assert_eq!(structure.atom_count(), 4);
    assert_eq!(structure.chain_count(), 2);
    assert!(!structure.is_ragged_ensemble());
}

#[test]
fn ragged_models_each_receive_perceived_bonds() {
    // Models that differ only in B-factor stay dense, so the second model changes
    // an atom's identity to force a different topology.
    let second = model(2, " 20.00").replace(" CA  GLY", " CB  GLY");
    let structure = read(format!("{}{second}END\n", model(1, " 10.00")));
    assert!(structure.is_ragged_ensemble());
    let Some(models) = structure.ragged_models() else {
        panic!("a ragged ensemble exposes its models")
    };
    assert_eq!(models.len(), 2);
    for model in models {
        assert_eq!(model.atom_count(), 4);
        assert_eq!(model.bonds().len(), 2, "one bond per chain");
    }
}

#[test]
fn a_real_nmr_ensemble_keeps_every_model_and_atom() {
    let Some(bytes) = molframe_bench::Sample::Ensemble.pdb() else {
        panic!("the ensemble sample ships a PDB file")
    };
    let Ok((structure, _)) = read_bytes(bytes.to_vec(), Some("2m7c.pdb"), &ReadOptions::new())
    else {
        panic!("the NMR ensemble must read")
    };
    assert_eq!(structure.model_count(), 32);
    assert_eq!(structure.atom_count(), 300);
}
