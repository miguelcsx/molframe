//! Scores for how alike two structures are.
//!
//! These take two coordinate sets that are already in correspondence — the
//! `i`-th point of one answers to the `i`-th of the other — and return a single
//! number. Deciding which atom matches which is a separate problem; these
//! measure a correspondence, they do not find one.
//!
//! The scores split into two kinds. lDDT is superposition-free: it compares the
//! internal distances of each structure, so no fitting is needed and no fitting
//! can bias it. TM-score and GDT first fit one structure onto the other with the
//! shared superposition kernel, then score the residuals.

#![forbid(unsafe_code)]

#[path = "equivalent.rs"]
mod atom_equivalence;
#[path = "ce.rs"]
mod combinatorial_extension;
#[path = "cad.rs"]
mod contact_area_difference;
#[path = "similarity.rs"]
mod contact_overlap;
#[path = "mapping/mod.rs"]
mod correspondence;
#[path = "dockq.rs"]
mod docking_quality;
#[path = "error.rs"]
mod failure;
mod governed_parameters;
mod interface;
#[path = "lddt.rs"]
mod local_distance;
mod numeric;
#[path = "governed/mod.rs"]
mod policy_execution;
#[path = "qs.rs"]
mod quaternary;
pub mod region;
pub mod superposed;
pub mod workflow;

pub use atom_equivalence::{
    EquivalentAtomMapping, LigandRmsd, equivalent_atom_mappings, ligand_symmetry_rmsd,
};
pub use combinatorial_extension::{
    CeAlignment, CeError, CeOptions, CeSignificanceProfile, ce_align, ce_alignments,
};
pub use contact_area_difference::{
    CadConstructionError, CadContact, CadError, CadScore, ContactArea, LocalCad, cad_contact_areas,
    cad_score,
};
pub use contact_overlap::{ContactSimilarity, contact_map_similarity};
pub use correspondence::{
    ChainAlternative, ChainAssignment, ChainMapping, ChainSequence, ResidueMatch, assign_chains,
    chain_sequences, map_chains, map_sequence_to_structure,
};
pub use docking_quality::{DockQ, DockQOptions, dockq, dockq_in_namespace};
pub use failure::CompareError;
pub use local_distance::{EmptyLddtPolicy, LddtOptions, lddt, lddt_with_options};
pub use policy_execution::{
    GovernedCompareError, governed_assign_chains, governed_cad_contact_areas, governed_cad_score,
    governed_ce_align, governed_ce_alignments, governed_comparison_workflow,
    governed_contact_map_similarity, governed_dockq, governed_equivalent_atom_mappings,
    governed_gdt_ha, governed_gdt_ts, governed_gdt_with_cutoffs, governed_interface_rmsd,
    governed_lddt, governed_ligand_symmetry_rmsd, governed_map_chains,
    governed_map_sequence_to_structure, governed_pocket_rmsd, governed_qs_score, governed_tm_score,
    governed_weighted_rmsd,
};
pub use quaternary::{EmptyQsPolicy, QsOptions, qs_score, qs_score_in_namespace};
pub use region::{InterfaceRmsd, PocketRmsd, interface_rmsd, pocket_rmsd};
pub use superposed::{gdt_ha, gdt_ts, gdt_with_cutoffs, tm_score, weighted_rmsd};
pub use workflow::{
    ComparisonAlignment, ComparisonVerdict, DistanceMeasurement, PointMapping, PointMatch,
    align_mapping, decide_rmsd, measure_mapping,
};
