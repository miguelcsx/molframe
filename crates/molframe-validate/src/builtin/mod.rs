//! Built-in validation against caller-selected reference tables.
//!
//! Everything here works from standard residue and atom names and the
//! reference tables a [`ReferenceSet`] resolves (a path, a cache entry or the
//! profile's choice), so no external library or annotation is needed to validate a
//! protein or nucleic acid structure. No table is compiled into the library.

#[path = "clashscore.rs"]
mod clash_scoring;
#[path = "error.rs"]
mod failure;
#[cfg(test)]
#[path = "fixtures.rs"]
mod fixtures;
#[path = "geometry_z.rs"]
mod geometry_scoring;
#[path = "nucleic_policy.rs"]
mod nucleic_policy;
#[path = "nucleic_conformation.rs"]
mod nucleic_scoring;
#[path = "protein.rs"]
mod protein_conformation;
#[path = "data.rs"]
mod reference_data;
#[path = "reference.rs"]
mod reference_set;
#[path = "tables.rs"]
mod reference_tables;
#[path = "layout.rs"]
mod residue_layout;
#[path = "summary.rs"]
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
