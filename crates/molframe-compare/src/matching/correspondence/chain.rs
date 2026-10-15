//! CCD-backed extraction and comparison of polymer chain sequences.

use molframe_chem::{ComponentKind, ComponentProvider};
use molframe_core::contract::Namespace;
use molframe_core::structure::Structure;
use molframe_core::{Code, Diagnostic};
use molframe_seq::{AlignError, Scoring, global};

use crate::numeric::u64_to_f64;

/// A chain reduced to its one-letter sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainSequence {
    /// The chain identifier.
    pub label: String,
    /// The one-letter sequence of its standard residues.
    pub sequence: Vec<u8>,
}

/// Extracts the one-letter sequence of every chain with standard residues.
///
/// # Errors
///
/// Returns a CCD provider diagnostic or reports a missing component/polymer code.
pub fn chain_sequences(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
) -> Result<Vec<ChainSequence>, Diagnostic> {
    let mut chains = Vec::new();
    for chain in structure.data().chains() {
        let mut sequence = Vec::new();
        for residue in chain.residues() {
            if let Some(name) = residue.name()
                && let Some(code) = component_code(provider, name)?
            {
                sequence.push(code);
            }
        }
        if sequence.is_empty() {
            continue;
        }
        let label = match chain_identifier(chain, namespace)? {
            Some(label) => label.to_string(),
            None => continue,
        };
        chains.push(ChainSequence { label, sequence });
    }
    Ok(chains)
}

/// Pairs two named chains and measures how identical their sequences are.
///
/// This is the single-pair counterpart of [`super::assign_chains`] for callers
/// who decide the pairing themselves. The identity is the one the assignment
/// would have used for the same two chains.
///
/// # Errors
///
/// Returns a diagnostic when a CCD component is missing, the namespace is not
/// supported, or either chain is not a polymer chain of the structure.
pub fn pair_chains(
    reference: &Structure,
    target: &Structure,
    labels: (&str, &str),
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    scoring: Scoring,
) -> Result<super::ChainMapping, Diagnostic> {
    let find = |structure: &Structure, label: &str| -> Result<Vec<u8>, Diagnostic> {
        chain_sequences(structure, provider, namespace)?
            .into_iter()
            .find(|chain| chain.label == label)
            .map(|chain| chain.sequence)
            .ok_or_else(|| {
                Diagnostic::new(Code::E4003)
                    .with_context("required", "polymer chain with standard residues")
                    .with_context("chain", label)
            })
    };
    let left = find(reference, labels.0)?;
    let right = find(target, labels.1)?;
    Ok(super::ChainMapping {
        reference: labels.0.to_owned(),
        target: labels.1.to_owned(),
        identity: identity(&left, &right, scoring)?,
    })
}

pub(super) fn chain_identifier(
    chain: molframe_core::structure::ChainRef<'_>,
    namespace: Namespace,
) -> Result<Option<&str>, Diagnostic> {
    match namespace {
        Namespace::Label => Ok(chain.label()),
        Namespace::Auth => Ok(chain.auth_label()),
        Namespace::Explicit => Err(Diagnostic::new(Code::E4003)
            .with_context("namespace", "explicit chain namespace required")),
        _ => Err(Diagnostic::new(Code::E4003).with_context("namespace", "unsupported")),
    }
}

pub(super) fn component_code(
    provider: &dyn ComponentProvider,
    component_id: &str,
) -> Result<Option<u8>, Diagnostic> {
    let component = provider.get(component_id)?.ok_or_else(|| {
        Diagnostic::new(Code::E4003)
            .with_context("required", "CCD component for sequence mapping")
            .with_context("component", component_id)
    })?;
    if !matches!(
        component.kind,
        ComponentKind::AminoAcid | ComponentKind::Nucleotide
    ) {
        return Ok(None);
    }
    component.one_letter_code.map(Some).ok_or_else(|| {
        Diagnostic::new(Code::E4003)
            .with_context("required", "CCD one-letter polymer code")
            .with_context("component", component_id)
    })
}

pub(super) fn identity(left: &[u8], right: &[u8], scoring: Scoring) -> Result<f64, Diagnostic> {
    let alignment = global(left, right, scoring).map_err(alignment_diagnostic)?;
    let mut aligned = 0u64;
    let mut matches = 0u64;
    for column in &alignment.columns {
        let (Some(left_index), Some(right_index)) = (column.left, column.right) else {
            continue;
        };
        aligned += 1;
        if left[left_index] == right[right_index] {
            matches += 1;
        }
    }
    if aligned == 0 {
        Ok(0.0)
    } else {
        Ok(u64_to_f64(matches) / u64_to_f64(aligned))
    }
}

pub(super) fn alignment_diagnostic(error: AlignError) -> Diagnostic {
    Diagnostic::new(Code::E4003).with_context("alignment", error.to_string())
}
