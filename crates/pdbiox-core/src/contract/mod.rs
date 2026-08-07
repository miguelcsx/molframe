//! What an analysis assumed, and what it is therefore worth.
//!
//! Every analysis takes a policy and returns a value plus its coverage, its
//! status, its warnings and its provenance. Decisions that are otherwise silent
//! defaults — which assembly, which model, which conformation, which identifier
//! namespace, what to do about atoms that were never modelled — become data that
//! can be inspected, fingerprinted, varied and published.

mod analysis;
mod policy;
mod provenance;
mod serialise;

pub use analysis::{Analysis, Assumption, AssumptionSource, Coverage, ImpactEstimate, Status};
pub use policy::{
    AlignmentPolicy, AltlocPolicy, AnalysisPolicy, AssemblyChoice, ContactDefinition,
    EquivalencePolicy, Fingerprint, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace,
    PeriodicPolicy, PolicyField, Precision, ProfileId, RadiiSet, SymmetryPolicy, Tolerance,
};
pub use provenance::{DictionaryVersion, Provenance, SourceRef};
