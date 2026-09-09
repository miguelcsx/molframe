use super::*;
use crate::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_atoms\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\nATOM 2 O O1 LIG A 1 1.3 0 0\n\
ATOM 3 C C2 LIG A 1 5 0 0\n";

#[test]
fn inference_adds_only_radius_compatible_nearby_pairs_with_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = pdbiox_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let report = infer_bonds(
        &structure,
        BondInference::default(),
        &pdbiox_core::ExecutionContext::default(),
    )
    .expect("inference works");
    let bonds: Vec<_> = report.structure.data().bonds.iter().collect();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].provenance, BondProvenance::InferredDistance);
    assert!(report.skipped_atoms.is_empty());
}
