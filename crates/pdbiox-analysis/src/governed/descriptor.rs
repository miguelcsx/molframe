//! Stable algorithm identity and sorted result-affecting parameters.

use pdbiox_core::contract::{AlgorithmId, AnalysisParameters, ParameterValue, Provenance};

/// Everything an adapter must add to algorithm provenance.
#[derive(Clone, Debug)]
pub struct AnalysisDescriptor {
    algorithm_name: Box<str>,
    algorithm_version: Box<str>,
    parameters: AnalysisParameters,
}

impl AnalysisDescriptor {
    /// Starts an explicitly versioned descriptor.
    #[must_use]
    pub fn new(name: impl Into<Box<str>>, version: impl Into<Box<str>>) -> Self {
        Self {
            algorithm_name: name.into(),
            algorithm_version: version.into(),
            parameters: AnalysisParameters::new(),
        }
    }

    /// Records or replaces one typed result-affecting parameter.
    #[must_use]
    pub fn with_parameter(mut self, name: impl Into<Box<str>>, value: ParameterValue) -> Self {
        self.parameters.insert(name.into(), value);
        self
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
        provenance
    }
}
