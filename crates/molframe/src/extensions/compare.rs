//! Structure-first comparison operations, as methods.
//!
//! Every method here forwards to a function in [`crate::compare`]. The receiver
//! is the structure the free function takes first: `model.dockq(native, …)` and
//! `reference.map_chains(target, …)`. Those two orderings are the kernels' own
//! and are kept rather than smoothed over, because a method that silently
//! swapped them would be the one thing a reader could not check.
//!
//! The coordinate-first comparators — lDDT, TM-score, GDT, the region RMSDs,
//! CAD and CE — are not here. They take `&[[f32; 3]]` and a mapping, so the
//! structure-level call is a gather the caller writes down, not a forward.

use crate::structure::Structure;
use crate::{AnalysisPolicy, Namespace};
use molframe_chem::ComponentProvider;
use molframe_seq::Scoring;

/// Comparison operations over a structure.
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "compare")]
/// # {
/// use molframe::prelude::*;
/// use molframe::compare::QsOptions;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (model, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let (native, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let score = model.qs_score(&native, "A", "A", QsOptions::standard(5.0));
/// assert!(score.is_ok());
/// # }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait CompareExt {
    /// The `DockQ` score of this model against `native`, over two chains.
    ///
    /// # Errors
    ///
    /// Returns a comparison error when a named chain is absent or the two
    /// structures share too few atoms to align.
    fn dockq(
        &self,
        native: &Structure,
        receptor: &str,
        ligand: &str,
        options: molframe_compare::DockQOptions,
    ) -> Result<molframe_compare::DockQ, molframe_compare::CompareError>;

    /// [`Self::dockq`] reading chain names in `namespace`.
    ///
    /// # Errors
    ///
    /// As [`Self::dockq`].
    fn dockq_in_namespace(
        &self,
        native: &Structure,
        receptor: &str,
        ligand: &str,
        namespace: Namespace,
        options: molframe_compare::DockQOptions,
    ) -> Result<molframe_compare::DockQ, molframe_compare::CompareError>;

    /// The QS score of this model against `native`, over two chains.
    ///
    /// # Errors
    ///
    /// As [`Self::dockq`].
    fn qs_score(
        &self,
        native: &Structure,
        first_chain: &str,
        second_chain: &str,
        options: molframe_compare::QsOptions,
    ) -> Result<f64, molframe_compare::CompareError>;

    /// [`Self::qs_score`] reading chain names in `namespace`.
    ///
    /// # Errors
    ///
    /// As [`Self::dockq`].
    fn qs_score_in_namespace(
        &self,
        native: &Structure,
        first_chain: &str,
        second_chain: &str,
        namespace: Namespace,
        options: molframe_compare::QsOptions,
    ) -> Result<f64, molframe_compare::CompareError>;

    /// Maps this reference's chains onto `target`'s.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the provider lacks a sequence for a chain.
    fn map_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
        min_identity: f64,
    ) -> Result<Vec<molframe_compare::ChainMapping>, molframe_core::diagnostic::Diagnostic>;

    /// Assigns this reference's chains onto `target`'s, one to one.
    ///
    /// # Errors
    ///
    /// As [`Self::map_chains`].
    fn assign_chains(
        &self,
        target: &Structure,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
        min_identity: f64,
    ) -> Result<molframe_compare::ChainAssignment, molframe_core::diagnostic::Diagnostic>;

    /// The sequence of every chain.
    ///
    /// # Errors
    ///
    /// As [`Self::map_chains`].
    fn chain_sequences(
        &self,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
    ) -> Result<Vec<molframe_compare::ChainSequence>, molframe_core::diagnostic::Diagnostic>;

    /// Aligns `query` against this structure's chains.
    ///
    /// # Errors
    ///
    /// As [`Self::map_chains`].
    fn map_sequence_to_structure(
        &self,
        query: &[u8],
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        scoring: Scoring,
    ) -> Result<Vec<molframe_compare::ResidueMatch>, molframe_core::diagnostic::Diagnostic>;

    /// [`Self::dockq`] carrying its policy and provenance.
    ///
    /// # Errors
    ///
    /// Returns a governed error for a missing chain or an unmet assumption.
    fn governed_dockq(
        &self,
        reference: &Structure,
        receptor: &str,
        ligand: &str,
        options: molframe_compare::DockQOptions,
        policy: &AnalysisPolicy,
    ) -> Result<crate::Analysis<molframe_compare::DockQ>, molframe_compare::GovernedCompareError>;

    /// [`Self::qs_score`] carrying its policy and provenance.
    ///
    /// # Errors
    ///
    /// As [`Self::governed_dockq`].
    fn governed_qs_score(
        &self,
        reference: &Structure,
        first_chain: &str,
        second_chain: &str,
        options: molframe_compare::QsOptions,
        policy: &AnalysisPolicy,
    ) -> Result<crate::Analysis<f64>, molframe_compare::GovernedCompareError>;

    /// [`Self::map_chains`] carrying its policy and provenance.
    ///
    /// # Errors
    ///
    /// As [`Self::governed_dockq`].
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
    >;

    /// [`Self::assign_chains`] carrying its policy and provenance.
    ///
    /// # Errors
    ///
    /// As [`Self::governed_dockq`].
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
    >;

    /// [`Self::map_sequence_to_structure`] carrying its policy and provenance.
    ///
    /// # Errors
    ///
    /// As [`Self::governed_dockq`].
    fn governed_map_sequence_to_structure(
        &self,
        query: &[u8],
        provider: &dyn ComponentProvider,
        scoring: Scoring,
        policy: &AnalysisPolicy,
    ) -> Result<
        crate::Analysis<Vec<molframe_compare::ResidueMatch>>,
        molframe_compare::GovernedCompareError,
    >;
}

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
