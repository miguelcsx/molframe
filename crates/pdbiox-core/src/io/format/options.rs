//! Shared parsing options and reader/writer traits.

use super::super::input::{InputBuffer, Limits};
use super::detect::Format;
use crate::diagnostic::Diagnostic;
use crate::diagnostic::{Diagnostics, Strictness};
use crate::structure::Structure;

/// How much irregularity a read tolerates, and what it does about it.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
#[non_exhaustive]
pub enum ParseMode {
    /// Anything violating the specification fails the read.
    Strict,
    /// Produce a structure and report what was wrong with it.
    #[default]
    Permissive,
    /// Continue past local errors, marking the regions they affected.
    Recover,
}

impl ParseMode {
    /// The strictness this mode implies.
    #[must_use]
    pub const fn strictness(self) -> Strictness {
        match self {
            Self::Strict => Strictness::Strict,
            Self::Permissive => Strictness::Medium,
            Self::Recover => Strictness::Loose,
        }
    }
}

/// Policy for atom rows whose format-specific element field is absent or invalid.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum MissingElementPolicy {
    /// Preserve the missing value as [`Element::UNKNOWN`].
    #[default]
    PreserveUnknown,
    /// Explicitly apply the naming convention of the selected input format.
    InferFromAtomName,
}

/// Policy for adjacent mmCIF atom groups whose residue identifiers are identical.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum AmbiguousResidueBoundaryPolicy {
    /// Refuse to invent a residue boundary that the deposited identifiers do not define.
    #[default]
    Reject,
    /// Explicitly split when an atom name repeats in the same alternate location.
    InferFromFileOrder,
}

/// What a read should and should not bother doing.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    /// Which format, or automatic detection.
    pub format: Format,
    /// How much irregularity to tolerate.
    pub mode: ParseMode,
    /// Read only the first model.
    pub only_first_model: bool,
    /// Read coordinates and skip everything else.
    pub only_atomic_coords: bool,
    /// Drop hydrogens while reading.
    pub discard_hydrogens: bool,
    /// What to do when a row does not declare a valid element.
    pub missing_element_policy: MissingElementPolicy,
    /// What to do when mmCIF residue identifiers do not define a boundary.
    pub ambiguous_residue_boundary_policy: AmbiguousResidueBoundaryPolicy,
    /// Ceilings that apply while reading.
    pub limits: Limits,
}

impl ReadOptions {
    /// Options that read everything, tolerating what the archive contains.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the format.
    #[must_use]
    pub const fn format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Sets how much irregularity to tolerate.
    #[must_use]
    pub const fn mode(mut self, mode: ParseMode) -> Self {
        self.mode = mode;
        self
    }

    /// Reads only the first model.
    #[must_use]
    pub const fn only_first_model(mut self, only: bool) -> Self {
        self.only_first_model = only;
        self
    }

    /// Reads coordinates and skips everything else.
    #[must_use]
    pub const fn only_atomic_coords(mut self, only: bool) -> Self {
        self.only_atomic_coords = only;
        self
    }

    /// Selects an explicit policy for missing or invalid element fields.
    #[must_use]
    pub const fn missing_element_policy(mut self, policy: MissingElementPolicy) -> Self {
        self.missing_element_policy = policy;
        self
    }

    /// Selects an explicit policy for ambiguous mmCIF residue boundaries.
    #[must_use]
    pub const fn ambiguous_residue_boundary_policy(
        mut self,
        policy: AmbiguousResidueBoundaryPolicy,
    ) -> Self {
        self.ambiguous_residue_boundary_policy = policy;
        self
    }

    /// Drops hydrogens while reading.
    #[must_use]
    pub const fn discard_hydrogens(mut self, discard: bool) -> Self {
        self.discard_hydrogens = discard;
        self
    }

    /// Sets the ceilings that apply while reading.
    #[must_use]
    pub const fn limits(mut self, limits: Limits) -> Self {
        self.limits = limits;
        self
    }

    /// Finishes a read under this option set.
    ///
    /// # Errors
    ///
    /// Returns the ordered findings when at least one finding is an error under
    /// the selected parse mode.
    pub fn finish(
        &self,
        structure: Structure,
        findings: impl IntoIterator<Item = Diagnostic>,
    ) -> ReadResult {
        let mut ordered = Diagnostics::new();
        ordered.extend(findings);
        let is_error = ordered.has_error(self.mode.strictness());
        let findings = ordered.finish();
        if is_error {
            Err(findings)
        } else {
            Ok((structure, findings))
        }
    }
}

/// A structure and everything that was wrong with the file it came from.
pub type ReadResult = Result<(Structure, Vec<Diagnostic>), Vec<Diagnostic>>;

/// A reader for one format.
pub trait Reader {
    /// The format this reads.
    const FORMAT: Format;

    /// Reads a structure.
    ///
    /// # Errors
    ///
    /// Returns the findings that stopped the read.
    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult;
}

/// Which parts of a structure a write should include.
pub trait Select {
    /// Whether to write this model.
    fn accept_model(&self, _model: crate::index::ModelIndex) -> bool {
        true
    }

    /// Whether to write this chain.
    fn accept_chain(&self, _chain: crate::index::ChainIndex) -> bool {
        true
    }

    /// Whether to write this residue.
    fn accept_residue(&self, _residue: crate::index::ResidueIndex) -> bool {
        true
    }

    /// Whether to write this atom.
    fn accept_atom(&self, _atom: crate::index::AtomIndex) -> bool {
        true
    }
}

/// A filter that writes everything.
#[derive(Clone, Copy, Debug, Default)]
pub struct SelectAll;

impl Select for SelectAll {}
