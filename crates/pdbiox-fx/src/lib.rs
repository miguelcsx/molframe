//! Declarative, chemistry-aware functional geometry.

#![forbid(unsafe_code)]

mod alignment;
mod benchmark;
mod config;
mod evaluation;
mod mapping;
mod measurement;
mod numeric;
mod profile;
mod specification;
mod verdict;

pub use alignment::{AlignedMotif, AlignmentKind, align_intrinsic, align_with_transform};
pub use benchmark::{
    AmeMeasurementError, AmeMeasurements, MotifBenchMeasurements, measure_ame, measure_motifbench,
};
pub use config::{EvaluationSpecification, SpecificationError, read_evaluation_specification};
pub use evaluation::{Evaluation, EvaluationError, EvaluationReport, evaluate_motif};
pub use mapping::{MappedMotif, MappingError, MappingSet, map_motif};
pub use measurement::{
    ConstraintMeasurement, IndeterminateReason, MeasurementError, MeasurementOptions,
    MeasurementSet, MeasurementValue, measure_constraints,
};
pub use profile::{
    CandidateAggregation, CompatibilityProfile, CompatibilityVerdict, ProfileAlignment,
    ProfileAtomSet, ProfileMetric, ame_heavy_atom_1_0, motifbench_1_0,
};
pub use specification::{
    AtomSite, ComponentRole, ComponentSpec, Constraint, Motif, MotifError, NamedConstraint,
};
pub use verdict::{
    Comparison, MissingVerdict, RuleOutcome, Verdict, VerdictProfile, VerdictRule, VerdictStatus,
};
