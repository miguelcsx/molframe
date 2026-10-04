//! Stable algorithm identity and sorted result-affecting parameters.

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
        if self.reads.is_empty() {
            provenance
        } else {
            provenance.with_policy_reads(&self.reads)
        }
    }
}
