use super::PyStructure;
use molframe::Structure;

#[test]
fn structure_binding_forwards_counts_without_reinterpreting_them() {
    let inner = Structure::new(molframe::StructureData::empty());
    let bound = PyStructure::new(inner.clone());
    assert_eq!(bound.atom_count(), inner.atom_count());
    assert_eq!(bound.model_count(), inner.model_count());
    assert_eq!(bound.chain_count(), inner.chain_count());
    assert_eq!(bound.residue_count(), inner.residue_count());
}
