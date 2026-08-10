//! Alignment of an explicit query sequence to structure residues.

use super::chain::{alignment_diagnostic, component_code};
use pdbiox_chem::ComponentProvider;
use pdbiox_core::contract::Namespace;
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::Structure;
use pdbiox_core::{Code, Diagnostic};
use pdbiox_seq::{Scoring, semi_global};

/// A query sequence position aligned to a structure residue.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResidueMatch {
    /// Position in the query sequence.
    pub query_position: usize,
    /// The structure residue it aligned to.
    pub residue: ResidueIndex,
}

/// Maps a query sequence onto a structure's residues by alignment.
///
/// # Errors
///
/// Returns a CCD provider diagnostic or reports a missing component/polymer code.
pub fn map_sequence_to_structure(
    query: &[u8],
    structure: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    scoring: Scoring,
) -> Result<Vec<ResidueMatch>, Diagnostic> {
    let mut residues = Vec::new();
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            if let Some(name) = residue_identifier(residue, namespace)?
                && let Some(code) = component_code(provider, name)?
            {
                residues.push((residue.index(), code));
            }
        }
    }
    let structure_sequence: Vec<u8> = residues.iter().map(|(_, code)| *code).collect();
    let alignment =
        semi_global(query, &structure_sequence, scoring).map_err(alignment_diagnostic)?;
    Ok(alignment
        .columns
        .iter()
        .filter_map(|column| {
            let (Some(query_position), Some(structure_index)) = (column.left, column.right) else {
                return None;
            };
            Some(ResidueMatch {
                query_position,
                residue: residues[structure_index].0,
            })
        })
        .collect())
}

fn residue_identifier(
    residue: pdbiox_core::structure::ResidueRef<'_>,
    namespace: Namespace,
) -> Result<Option<&str>, Diagnostic> {
    match namespace {
        Namespace::Label => Ok(residue.name()),
        Namespace::Auth => Ok(residue.auth_name()),
        Namespace::Explicit => Err(Diagnostic::new(Code::E4003)
            .with_context("namespace", "explicit residue namespace required")),
        _ => Err(Diagnostic::new(Code::E4003).with_context("namespace", "unsupported")),
    }
}
