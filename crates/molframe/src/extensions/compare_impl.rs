use super::compare::CompareExt;
use crate::structure::Structure;
use crate::{AnalysisPolicy, Namespace};
use molframe_chem::ComponentProvider;
use molframe_seq::Scoring;

impl CompareExt for Structure {
    fn dockq(
        &self,
        native: &Structure,
        receptor: &str,
        ligand: &str,
        options: molframe_compare::DockQOptions,
    ) -> Result<molframe_compare::DockQ, molframe_compare::CompareError> {
        molframe_compare::dockq(self.engine(), native.engine(), receptor, ligand, options)
    }

    fn dockq_in_namespace(
        &self,
        native: &Structure,
        receptor: &str,
        ligand: &str,
        namespace: Namespace,
        options: molframe_compare::DockQOptions,
    ) -> Result<molframe_compare::DockQ, molframe_compare::CompareError> {
        molframe_compare::dockq_in_namespace(
            self.engine(),
            native.engine(),
            receptor,
            ligand,
            namespace,
            options,
        )
    }

    fn qs_score(
        &self,
        native: &Structure,
        first_chain: &str,
        second_chain: &str,
        options: molframe_compare::QsOptions,
    ) -> Result<f64, molframe_compare::CompareError> {
        molframe_compare::qs_score(
            self.engine(),
            native.engine(),
            first_chain,
            second_chain,
            options,
        )
    }

    fn qs_score_in_namespace(
        &self,
        native: &Structure,
        first_chain: &str,
        second_chain: &str,
        namespace: Namespace,
        options: molframe_compare::QsOptions,
    ) -> Result<f64, molframe_compare::CompareError> {
        molframe_compare::qs_score_in_namespace(
            self.engine(),
            native.engine(),
            first_chain,
            second_chain,
            namespace,
            options,
        )
    }

    fn mapped_dockq(
        &self,
        native: &Structure,
        receptor: &str,
        ligand: &str,
        mapping: &molframe_compare::MappingOptions<'_>,
        options: molframe_compare::DockQOptions,
    ) -> Result<
        (molframe_compare::DockQ, molframe_compare::MappedComparison),
        molframe_compare::MappedCompareError,
    > {
        molframe_compare::mapped_dockq(
            self.engine(),
            native.engine(),
            receptor,
            ligand,
            mapping,
            options,
        )
    }

    fn mapped_qs_score(
        &self,
        native: &Structure,
        first_chain: &str,
        second_chain: &str,
        mapping: &molframe_compare::MappingOptions<'_>,
        options: molframe_compare::QsOptions,
    ) -> Result<(f64, molframe_compare::MappedComparison), molframe_compare::MappedCompareError>
    {
        molframe_compare::mapped_qs_score(
            self.engine(),
            native.engine(),
            first_chain,
            second_chain,
            mapping,
            options,
        )
    }

    fn map_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
        min_identity: f64,
    ) -> Result<Vec<molframe_compare::ChainMapping>, molframe_core::diagnostic::Diagnostic> {
        molframe_compare::map_chains(
            self.engine(),
            target.engine(),
            provider,
            namespace,
            scoring,
            min_identity,
        )
    }

    fn assign_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
        min_identity: f64,
    ) -> Result<molframe_compare::ChainAssignment, molframe_core::diagnostic::Diagnostic> {
        molframe_compare::assign_chains(
            self.engine(),
            target.engine(),
            provider,
            namespace,
            scoring,
            min_identity,
        )
    }

    fn chain_sequences(
        &self,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
    ) -> Result<Vec<molframe_compare::ChainSequence>, molframe_core::diagnostic::Diagnostic> {
        molframe_compare::chain_sequences(self.engine(), provider, namespace)
    }

    fn map_sequence_to_structure(
        &self,
        query: &[u8],
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
    ) -> Result<Vec<molframe_compare::ResidueMatch>, molframe_core::diagnostic::Diagnostic> {
        molframe_compare::map_sequence_to_structure(
            query,
            self.engine(),
            provider,
            namespace,
            scoring,
        )
    }

    fn governed_dockq(
        &self,
        reference: &Structure,
        receptor: &str,
        ligand: &str,
        options: molframe_compare::DockQOptions,
        policy: &AnalysisPolicy,
    ) -> Result<crate::Analysis<molframe_compare::DockQ>, molframe_compare::GovernedCompareError>
    {
        molframe_compare::governed_dockq(
            self.engine(),
            reference.engine(),
            receptor,
            ligand,
            options,
            policy,
        )
    }

    fn governed_qs_score(
        &self,
        reference: &Structure,
        first_chain: &str,
        second_chain: &str,
        options: molframe_compare::QsOptions,
        policy: &AnalysisPolicy,
    ) -> Result<crate::Analysis<f64>, molframe_compare::GovernedCompareError> {
        molframe_compare::governed_qs_score(
            self.engine(),
            reference.engine(),
            first_chain,
            second_chain,
            options,
            policy,
        )
    }

    fn governed_map_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        scoring: Scoring,
        minimum_identity: f64,
        policy: &AnalysisPolicy,
    ) -> Result<
        crate::Analysis<Vec<molframe_compare::ChainMapping>>,
        molframe_compare::GovernedCompareError,
    > {
        molframe_compare::governed_map_chains(
            self.engine(),
            target.engine(),
            provider,
            scoring,
            minimum_identity,
            policy,
        )
    }

    fn governed_assign_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        scoring: Scoring,
        minimum_identity: f64,
        policy: &AnalysisPolicy,
    ) -> Result<
        crate::Analysis<molframe_compare::ChainAssignment>,
        molframe_compare::GovernedCompareError,
    > {
        molframe_compare::governed_assign_chains(
            self.engine(),
            target.engine(),
            provider,
            scoring,
            minimum_identity,
            policy,
        )
    }

    fn governed_map_sequence_to_structure(
        &self,
        query: &[u8],
        provider: &dyn ComponentProvider,
        scoring: Scoring,
        policy: &AnalysisPolicy,
    ) -> Result<
        crate::Analysis<Vec<molframe_compare::ResidueMatch>>,
        molframe_compare::GovernedCompareError,
    > {
        molframe_compare::governed_map_sequence_to_structure(
            query,
            self.engine(),
            provider,
            scoring,
            policy,
        )
    }
}
