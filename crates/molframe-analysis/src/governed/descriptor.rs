//! Stable algorithm identity and sorted result-affecting parameters.

use super::{ForbiddenResolution, Requirement};
use molframe_core::contract::{
    AlgorithmId, AnalysisParameters, ParameterValue, PolicyField, Provenance,
};

/// Everything an adapter must add to algorithm provenance.
#[derive(Clone, Debug)]
pub struct AnalysisDescriptor {
    algorithm_name: Box<str>,
    algorithm_version: Box<str>,
    parameters: AnalysisParameters,
    reads: Vec<PolicyField>,
    replicated_systems: bool,
    forbidden: Vec<ForbiddenResolution>,
    required: Vec<Requirement>,
    estimand: Option<&'static str>,
}

impl AnalysisDescriptor {
    /// Starts an explicitly versioned descriptor.
    #[must_use]
    pub fn new(name: impl Into<Box<str>>, version: impl Into<Box<str>>) -> Self {
        Self {
            algorithm_name: name.into(),
            algorithm_version: version.into(),
            parameters: AnalysisParameters::new(),
            reads: Vec::new(),
            replicated_systems: true,
            forbidden: Vec::new(),
            required: Vec::new(),
            estimand: None,
        }
    }

    /// Records or replaces one typed result-affecting parameter.
    #[must_use]
    pub fn with_parameter(mut self, name: impl Into<Box<str>>, value: ParameterValue) -> Self {
        self.parameters.insert(name.into(), value);
        self
    }

    /// Declares policy fields this analysis applies beyond the ones the executor
    /// applies for every kernel.
    ///
    /// A field not declared here and not applied by the executor cannot change
    /// the result, and the audit refuses to vary it.
    #[must_use]
    pub fn reading(mut self, fields: &[PolicyField]) -> Self {
        self.reads.extend_from_slice(fields);
        self
    }

    /// Declares that atom correspondence across copies is ill-defined for this
    /// analysis, so it refuses to run over an assembly or crystal contacts.
    #[must_use]
    pub const fn without_replicated_systems(mut self) -> Self {
        self.replicated_systems = false;
        self
    }

    /// Declares a resolution of a decision under which this analysis has nothing to
    /// measure, so the executor refuses it rather than return an empty answer.
    #[must_use]
    pub fn forbidding(mut self, resolution: ForbiddenResolution) -> Self {
        self.forbidden.push(resolution);
        self
    }

    /// Declares information the analysed atoms must carry for there to be an answer;
    /// an input without it is indeterminate, not empty.
    #[must_use]
    pub fn requiring(mut self, requirement: Requirement) -> Self {
        self.required.push(requirement);
        self
    }

    /// States, in words, the quantity this analysis estimates.
    #[must_use]
    pub const fn estimating(mut self, estimand: &'static str) -> Self {
        self.estimand = Some(estimand);
        self
    }

    /// The resolutions this analysis refuses.
    #[must_use]
    pub fn forbidden_resolutions(&self) -> &[ForbiddenResolution] {
        &self.forbidden
    }

    /// The information this analysis needs.
    #[must_use]
    pub fn required_information(&self) -> &[Requirement] {
        &self.required
    }

    /// The quantity this analysis estimates, when it says.
    #[must_use]
    pub const fn estimand(&self) -> Option<&'static str> {
        self.estimand
    }

    /// The policy fields the descriptor declares, beyond the executor's own.
    #[must_use]
    pub fn policy_reads(&self) -> &[PolicyField] {
        &self.reads
    }

    /// Whether the analysis may run over a system with replicated atoms.
    #[must_use]
    pub const fn allows_replicated_systems(&self) -> bool {
        self.replicated_systems
    }

    /// Stable algorithm name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.algorithm_name
    }

    pub(super) fn apply(&self, mut provenance: Provenance) -> Provenance {
        provenance = provenance.with_algorithm(AlgorithmId::new(
            self.algorithm_name.clone(),
            self.algorithm_version.clone(),
        ));
        for (name, value) in &self.parameters {
            provenance = provenance.with_parameter(name.clone(), value.clone());
        }
        if let Some(estimand) = self.estimand {
            provenance = provenance.with_estimand(estimand);
        }
        if self.reads.is_empty() {
            provenance
        } else {
            provenance.with_policy_reads(&self.reads)
        }
    }
}
