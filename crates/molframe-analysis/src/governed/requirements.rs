//! What an analysis needs in order to have an answer, and what it cannot survive.
//!
//! A default that resolves a decision first can destroy the information an analysis
//! measures: a hydrogen-bond geometry has no answer once the hydrogens are gone, and
//! an empty table is a different statement from "there are none". Each analysis says
//! which resolutions it forbids and which information it requires.

use molframe_core::AtomSelection;
use molframe_core::contract::{AnalysisPolicy, PolicyField};
use molframe_core::structure::{AtomRef, Structure};

/// A resolution of one decision under which an analysis has nothing to measure.
#[derive(Clone, Copy, Debug)]
pub struct ForbiddenResolution {
    pub(super) field: PolicyField,
    pub(super) reason: &'static str,
    pub(super) applies: fn(&AnalysisPolicy) -> bool,
}

impl ForbiddenResolution {
    /// Forbids the values of `field` for which `applies` holds, saying why.
    #[must_use]
    pub const fn new(
        field: PolicyField,
        reason: &'static str,
        applies: fn(&AnalysisPolicy) -> bool,
    ) -> Self {
        Self {
            field,
            reason,
            applies,
        }
    }

    /// The decision whose value is refused.
    #[must_use]
    pub const fn field(&self) -> PolicyField {
        self.field
    }

    /// Why that value leaves nothing to measure.
    #[must_use]
    pub const fn reason(&self) -> &'static str {
        self.reason
    }

    /// Whether the policy chooses a value this analysis forbids.
    #[must_use]
    pub fn forbids(&self, policy: &AnalysisPolicy) -> bool {
        (self.applies)(policy)
    }
}

/// Information the input must carry for an analysis to have an answer.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Requirement {
    /// At least one hydrogen atom is modelled.
    ExplicitHydrogens,
}

impl Requirement {
    /// What is required, as the reason an input without it has no answer.
    #[must_use]
    pub const fn reason(self) -> &'static str {
        match self {
            Self::ExplicitHydrogens => {
                "the input models no hydrogen atoms, which this analysis measures"
            }
        }
    }

    /// Whether `structure` carries what is required.
    #[must_use]
    pub fn satisfied_by(self, structure: &Structure) -> bool {
        match self {
            Self::ExplicitHydrogens => structure.data().atoms().any(is_hydrogen),
        }
    }
}

fn is_hydrogen(atom: AtomRef<'_>) -> bool {
    atom.element()
        .is_some_and(|element| element.atomic_number() == 1)
}

/// The selection without its hydrogen atoms.
pub(super) fn without_hydrogens(structure: &Structure, chosen: &AtomSelection) -> AtomSelection {
    let mut hydrogen = vec![false; structure.atom_count() as usize];
    for atom in structure.data().atoms() {
        if let Some(slot) = hydrogen.get_mut(atom.index().as_usize()) {
            *slot = is_hydrogen(atom);
        }
    }
    AtomSelection::from_sorted(
        chosen
            .iter()
            .filter(|&atom| !matches!(hydrogen.get(atom as usize), Some(true)))
            .collect(),
    )
}
