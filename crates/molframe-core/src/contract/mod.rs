//! What an analysis assumed, and what it is therefore worth.
//!
//! Every analysis takes a policy and returns a value plus its coverage, its
//! status, its warnings and its provenance. Decisions that are otherwise silent
//! defaults — which assembly, which model, which conformation, which identifier
//! namespace, what to do about atoms that were never modelled — become data that
//! can be inspected, fingerprinted, varied and published.

mod analysis;
mod digest;
mod outcome;
mod policy;
mod provenance;
mod reexecution;
mod serialise;

pub use analysis::{
    Analysis, Assumption, AssumptionSource, Coverage, Impact, MeasuredImpact, Status,
};
pub use digest::ContentDigest;
pub use outcome::{Indeterminacy, MissingPolicyError, Outcome, Quality, resolve_missing};
pub(crate) use policy::closed_vocabulary;
pub use policy::{
    AlignmentPolicy, AltlocPolicy, AnalysisPolicy, AssemblyChoice, ContactDefinition,
    EquivalencePolicy, Fingerprint, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace,
    PeriodicPolicy, PolicyField, PolicyParseError, Precision, ProfileId, RadiiSet, SymmetryPolicy,
    Tolerance, canonical_spelling,
};
pub use provenance::{
    AlgorithmId, AnalysisParameters, DictionaryVersion, ParameterValue, Provenance, SourceRef,
};
pub use reexecution::{
    Reexecution, ReexecutionEnvironment, ReexecutionError, reexecute_from_provenance,
};
pub use serialise::push_json_string;
