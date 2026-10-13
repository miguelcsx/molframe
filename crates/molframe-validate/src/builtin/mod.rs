//! Built-in validation against bundled reference data.
//!
//! Everything here works from standard residue and atom names and the
//! reference tables under `data/`, so no external library or annotation is
//! needed to validate a protein or nucleic acid structure.

#[path = "bundled.rs"]
mod bundled_tables;
#[path = "clashscore.rs"]
mod clash_scoring;
#[path = "error.rs"]
mod failure;
#[path = "geometry_z.rs"]
mod geometry_scoring;
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
#[path = "user_tables.rs"]
mod user_reference;

pub use bundled_tables::bundled_reference;
pub use clash_scoring::{ClashScore, ClashScoreOptions, HydrogenMode, clash_score};
pub use failure::BuiltinError;
pub use geometry_scoring::{GeometryOutlier, GeometryScores, NO_ATOM, geometry_scores};
pub use protein_conformation::{NOT_EVALUATED, OmegaClass, ProteinTable, protein_table};
pub use reference_data::{
    RAMA_BINS, RAMA_FAVOURED, RamaCategory, RamaGrids, RamaThresholds, RotamerClass, RotamerTable,
    RotamerTables, TableError,
};
pub use reference_set::{Part, ReferenceOverrides, ReferenceSet, TableInfo};
pub use reference_tables::{ChiDefinition, ChiDefinitions, GeometryTargets, Target};
pub use residue_layout::Layout;
