//! What kind of uncertainty a decision carries.
//!
//! A choice between conformations and a choice of numerical tolerance are not the
//! same sort of doubt, and an audit that pools them reports a number that means
//! nothing to either audience. Each policy decision belongs to one class, and the
//! share of the variation a class explains can be read off separately.

use molframe_core::contract::PolicyField;
use std::fmt;

/// The kind of uncertainty a decision represents.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[non_exhaustive]
pub enum UncertaintyClass {
    /// Which structures and models exist to be analysed.
    Structural,
    /// How the deposited structure is read: which arrangement, which conformer,
    /// which identifiers, what a missing atom means, how hydrogens enter, which
    /// atoms correspond.
    Interpretive,
    /// Which method answers the question: how a contact is defined, which radii
    /// the distance is measured against, how structures are aligned.
    Algorithmic,
    /// Floating-point width and tolerance: the part of the answer that should not
    /// move and does.
    Numerical,
}

impl UncertaintyClass {
    /// Every class, in a fixed order.
    pub const ALL: [Self; 4] = [
        Self::Structural,
        Self::Interpretive,
        Self::Algorithmic,
        Self::Numerical,
    ];

    /// The class a policy decision belongs to.
    #[must_use]
    pub const fn of(field: PolicyField) -> Self {
        match field {
            PolicyField::Model | PolicyField::Periodic => Self::Structural,
            PolicyField::Assembly
            | PolicyField::Symmetry
            | PolicyField::Altloc
            | PolicyField::Identifiers
            | PolicyField::MissingAtoms
            | PolicyField::Hydrogens
            | PolicyField::AtomEquivalence => Self::Interpretive,
            PolicyField::Alignment | PolicyField::VdwRadii | PolicyField::ContactDef => {
                Self::Algorithmic
            }
            PolicyField::Precision | PolicyField::FloatTolerance => Self::Numerical,
            // A decision added later is a decision about how the answer is made until
            // someone classifies it.
            _ => Self::Algorithmic,
        }
    }

    /// The class's name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Structural => "structural",
            Self::Interpretive => "interpretive",
            Self::Algorithmic => "algorithmic",
            Self::Numerical => "numerical",
        }
    }
}

impl fmt::Display for UncertaintyClass {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}
