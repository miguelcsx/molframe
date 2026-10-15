//! Residue correspondence by identifier, for models that share a numbering.

use super::chain::chain_identifier;
use crate::ResidueCorrespondence;
use molframe_chem::{ComponentKind, ComponentProvider};
use molframe_core::Diagnostic;
use molframe_core::contract::Namespace;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::{ResidueRef, Structure};
use std::collections::BTreeMap;

/// Pairs every standard polymer residue of `reference` with the residue of
/// `model` that has the same chain, number, insertion code and component.
///
/// The identifier is the label chain and sequence position under
/// [`Namespace::Label`], and the author chain, residue number and insertion
/// code under [`Namespace::Auth`]. A reference residue the model lacks, or
/// holds under another component name, is paired with `None` so that the score
/// counts it as missing. Residues without a sequence identifier in the
/// namespace, and components the provider does not describe as amino acids or
/// nucleotides, are outside the domain. Reference order is kept.
///
/// # Errors
///
/// Returns a provider diagnostic or an unsupported-namespace diagnostic.
pub fn pair_residues_by_identifier(
    reference: &Structure,
    model: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
) -> Result<Vec<ResidueCorrespondence>, Diagnostic> {
    let mut available: BTreeMap<(&str, i32, &str), (ResidueIndex, &str)> = BTreeMap::new();
    for chain in model.data().chains() {
        let Some(label) = chain_identifier(chain, namespace)? else {
            continue;
        };
        for residue in chain.residues() {
            if let (Some((number, insertion)), Some(name)) =
                (identity(residue, namespace), residue.name())
            {
                available
                    .entry((label, number, insertion))
                    .or_insert((residue.index(), name));
            }
        }
    }

    let mut pairs = Vec::new();
    for chain in reference.data().chains() {
        let Some(label) = chain_identifier(chain, namespace)? else {
            continue;
        };
        for residue in chain.residues() {
            let (Some((number, insertion)), Some(name)) =
                (identity(residue, namespace), residue.name())
            else {
                continue;
            };
            if !is_polymer_residue(provider, name)? {
                continue;
            }
            let counterpart = match available.get(&(label, number, insertion)) {
                Some((index, found)) if *found == name => Some(*index),
                _ => None,
            };
            pairs.push(ResidueCorrespondence {
                reference: residue.index(),
                model: counterpart,
            });
        }
    }
    Ok(pairs)
}

fn identity(residue: ResidueRef<'_>, namespace: Namespace) -> Option<(i32, &str)> {
    match namespace {
        Namespace::Auth => Some((
            residue.auth_seq_id()?,
            match residue.ins_code() {
                Some(code) => code,
                None => "",
            },
        )),
        _ => Some((residue.label_seq_id()?, "")),
    }
}

fn is_polymer_residue(provider: &dyn ComponentProvider, name: &str) -> Result<bool, Diagnostic> {
    Ok(match provider.get(name)? {
        Some(component) => matches!(
            component.kind,
            ComponentKind::AminoAcid | ComponentKind::Nucleotide
        ),
        None => false,
    })
}
