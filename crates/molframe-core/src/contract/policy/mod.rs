//! What an analysis was told to assume.

mod fields;
mod fingerprint;
mod profile;
mod vocabulary;

pub use fields::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy, PolicyField, Precision,
    ProfileId, RadiiSet, SymmetryPolicy, Tolerance,
};
pub use fingerprint::Fingerprint;
pub use profile::AnalysisPolicy;
pub use vocabulary::{PolicyParseError, canonical_spelling};
