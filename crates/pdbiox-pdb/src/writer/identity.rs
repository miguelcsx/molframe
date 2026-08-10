//! Namespace-coherent PDB identifier projection.

use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef};

/// Identifier namespace projected into legacy PDB fields.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PdbIdentifierNamespace {
    /// Normalized `label_*` identifiers.
    #[default]
    Label,
    /// Depositor-provided `auth_*` identifiers.
    Auth,
}

impl PdbIdentifierNamespace {
    pub(super) const fn as_str(self) -> &'static str {
        match self {
            Self::Label => "label",
            Self::Auth => "auth",
        }
    }
}

pub(crate) fn chain_name(chain: ChainRef<'_>, namespace: PdbIdentifierNamespace) -> Option<&str> {
    match namespace {
        PdbIdentifierNamespace::Label => chain.label(),
        PdbIdentifierNamespace::Auth => chain.auth_label(),
    }
}

pub(crate) fn residue_sequence(
    residue: ResidueRef<'_>,
    namespace: PdbIdentifierNamespace,
) -> Option<i32> {
    match namespace {
        PdbIdentifierNamespace::Label => residue.label_seq_id(),
        PdbIdentifierNamespace::Auth => residue.auth_seq_id(),
    }
}

pub(crate) fn atom_name(atom: AtomRef<'_>, namespace: PdbIdentifierNamespace) -> Option<&str> {
    match namespace {
        PdbIdentifierNamespace::Label => atom.name(),
        PdbIdentifierNamespace::Auth => atom.auth_name(),
    }
}

pub(crate) fn component_name<'a>(
    atom: AtomRef<'a>,
    residue: ResidueRef<'a>,
    namespace: PdbIdentifierNamespace,
) -> Option<&'a str> {
    match namespace {
        PdbIdentifierNamespace::Label => atom.component_name(),
        PdbIdentifierNamespace::Auth => residue.auth_name(),
    }
}
