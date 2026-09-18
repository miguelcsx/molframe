//! Policy-bound comparison entry points with complete provenance.

mod common;
mod equivalence;
mod interfaces;
mod mapping;
mod scores;

pub use common::GovernedCompareError;
pub use equivalence::{governed_equivalent_atom_mappings, governed_ligand_symmetry_rmsd};
pub use interfaces::{governed_dockq, governed_qs_score};
pub use mapping::{
    governed_assign_chains, governed_comparison_workflow, governed_interface_rmsd,
    governed_map_chains, governed_map_sequence_to_structure, governed_pocket_rmsd,
};
pub use scores::{
    governed_cad_contact_areas, governed_cad_score, governed_ce_align, governed_ce_alignments,
    governed_contact_map_similarity, governed_gdt_ha, governed_gdt_ts, governed_gdt_with_cutoffs,
    governed_lddt, governed_tm_score, governed_weighted_rmsd,
};

#[cfg(test)]
#[path = "../governed_tests.rs"]
mod tests;
