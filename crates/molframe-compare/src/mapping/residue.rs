//! Residue correspondence inside each mapped chain pair.

use super::assignment::ChainMapping;
use super::chain::{alignment_diagnostic, chain_identifier, component_code};
use molframe_chem::ComponentProvider;
use molframe_core::Diagnostic;
use molframe_core::contract::Namespace;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::Structure;
use molframe_seq::{Scoring, global};

/// A reference residue and the target residue it aligns to.
pub type ResiduePair = (ResidueIndex, ResidueIndex);

/// Aligns the polymer residues of every mapped chain pair, residue by residue.
///
/// Each pair is globally aligned on the one-letter codes the provider assigns
/// to its standard residues, and only columns holding a residue on both sides
/// become a correspondence. Pairs are returned chain by chain in mapping
/// order, then in sequence order.
///
/// # Errors
///
/// Returns a provider, namespace or alignment diagnostic.
pub fn map_residues(
    reference: &Structure,
    target: &Structure,
    chains: &[ChainMapping],
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    scoring: Scoring,
) -> Result<Vec<ResiduePair>, Diagnostic> {
    let mut pairs = Vec::new();
    for mapping in chains {
        let left = chain_residues(reference, &mapping.reference, provider, namespace)?;
        let right = chain_residues(target, &mapping.target, provider, namespace)?;
        let left_codes: Vec<u8> = left.iter().map(|(_, code)| *code).collect();
        let right_codes: Vec<u8> = right.iter().map(|(_, code)| *code).collect();
        let alignment = global(&left_codes, &right_codes, scoring).map_err(alignment_diagnostic)?;
        for column in &alignment.columns {
            let (Some(first), Some(second)) = (column.left, column.right) else {
                continue;
            };
            if let (Some((reference, _)), Some((target, _))) = (left.get(first), right.get(second))
            {
                pairs.push((*reference, *target));
            }
        }
    }
    Ok(pairs)
}

/// The standard residues of the named chain with their one-letter codes.
fn chain_residues(
    structure: &Structure,
    label: &str,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
) -> Result<Vec<(ResidueIndex, u8)>, Diagnostic> {
    let mut residues = Vec::new();
    for chain in structure.data().chains() {
        if chain_identifier(chain, namespace)? != Some(label) {
            continue;
        }
        for residue in chain.residues() {
            if let Some(name) = residue.name()
                && let Some(code) = component_code(provider, name)?
            {
                residues.push((residue.index(), code));
            }
        }
    }
    Ok(residues)
}
