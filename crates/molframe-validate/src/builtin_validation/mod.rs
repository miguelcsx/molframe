//! Built-in validation against caller-selected reference tables.
//!
//! Everything here works from standard residue and atom names and the
//! reference tables a [`ReferenceSet`] resolves (a path, a cache entry or the
//! profile's choice), so no external library or annotation is needed to validate a
//! protein or nucleic acid structure. No table is compiled into the library.

mod clash_scoring;
mod failure;
#[cfg(test)]
mod fixtures;
mod geometry_scoring;
mod nucleic_policy;
mod nucleic_scoring;
mod protein_conformation;
mod reference_data;
mod reference_set;
mod reference_tables;
mod residue_layout;
mod summaries;

pub use clash_scoring::{ClashScore, ClashScoreOptions, HydrogenMode, clash_score};
pub use failure::BuiltinError;
pub use geometry_scoring::{GeometryOutlier, GeometryScores, NO_ATOM, geometry_scores};
pub use nucleic_policy::{NucleicPolicy, Window};
pub use nucleic_scoring::{
    NUCLEIC_TORSIONS, NucleicTable, PUCKER_C2_ENDO, PUCKER_C3_ENDO, PUCKER_OTHER, SUITE_DELTA,
    SUITE_GAMMA, SUITE_PREVIOUS_DELTA, nucleic_table,
};
pub use protein_conformation::{NOT_EVALUATED, OmegaClass, ProteinTable, protein_table};
pub use reference_data::{
    RAMA_BINS, RamaCategory, RamaGrids, RamaThresholds, RotamerClass, RotamerTable, RotamerTables,
};
pub use reference_set::{Part, ReferenceSet};
pub use reference_tables::{ChiDefinition, ChiDefinitions, GeometryTargets, Target};
pub use residue_layout::Layout;
pub use summaries::{NucleicSummary, ProteinSummary, percent};
