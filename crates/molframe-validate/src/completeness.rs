//! Completeness of each chain against its canonical sequence.
//!
//! A polymer's entity declares the sequence that should be there; the model
//! often has less — disordered loops, unresolved termini. This reports, per
//! chain, how much of the canonical sequence was observed and exactly which
//! canonical positions are missing, so a gap is named rather than merely counted.
//!
//! Chains with no canonical sequence — ligands, water, or a polymer whose entity
//! carries no sequence — are omitted, since there is nothing to be complete
//! against.

use molframe_core::contract::Namespace;
use molframe_core::structure::{ChainRef, ChainSequenceExt, MissingResidue, Structure};

/// How complete one chain is relative to its canonical sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainCompleteness {
    /// The chain identifier.
    pub chain: String,
    /// Number of residues actually modelled.
    pub observed: usize,
    /// Number of residues the canonical sequence expects.
    pub canonical: usize,
    /// The canonical positions with no observed residue.
    pub missing: Vec<MissingResidue>,
}

/// Why chain completeness could not identify its report domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CompletenessError {
    /// `Explicit` requires a caller to select one concrete identifier namespace.
    #[error("chain completeness requires label or auth identifier namespace")]
    ExplicitNamespace,
    /// A chain in the assessed polymer domain lacks its requested identifier.
    #[error("a polymer chain has no identifier in the requested namespace")]
    MissingChainIdentifier,
    /// This crate does not understand a newer namespace variant.
    #[error("unsupported chain identifier namespace")]
    UnsupportedNamespace,
}

/// Reports the completeness of every chain that has a canonical sequence.
///
/// Chains are returned in topology order.
///
/// Runs in `O(residues)` time.
///
/// # Errors
///
/// Returns an error for an explicit or unsupported namespace, or when a
/// polymer chain lacks its requested identifier.
pub fn completeness(
    structure: &Structure,
    namespace: Namespace,
) -> Result<Vec<ChainCompleteness>, CompletenessError> {
    let mut report = Vec::new();
    for chain in structure.data().chains() {
        let canonical = chain.canonical_sequence();
        if canonical.is_empty() {
            continue;
        }
        let label = chain_identifier(chain, namespace)?
            .ok_or(CompletenessError::MissingChainIdentifier)?
            .to_string();
        report.push(ChainCompleteness {
            chain: label,
            observed: chain.observed_sequence().len(),
            canonical: canonical.len(),
            missing: chain.missing_residues(),
        });
    }
    Ok(report)
}

fn chain_identifier(
    chain: ChainRef<'_>,
    namespace: Namespace,
) -> Result<Option<&str>, CompletenessError> {
    match namespace {
        Namespace::Label => Ok(chain.label()),
        Namespace::Auth => Ok(chain.auth_label()),
        Namespace::Explicit => Err(CompletenessError::ExplicitNamespace),
        _ => Err(CompletenessError::UnsupportedNamespace),
    }
}

#[cfg(test)]
#[path = "completeness_tests.rs"]
mod tests;
