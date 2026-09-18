//! Native deterministic SMARTS substructure queries over CCD components.

use crate::{Component, StereoConfiguration};
use molframe_core::{BondOrder, Element, Structure};
use std::fmt;

/// A parsed SMARTS query graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmartsPattern {
    pub(crate) atoms: Vec<AtomExpression>,
    pub(crate) bonds: Vec<PatternBond>,
}

impl SmartsPattern {
    /// Parses a SMARTS query without consulting a chemical dictionary.
    ///
    /// # Errors
    ///
    /// Returns [`SmartsError`] when the graph syntax is malformed or contains
    /// a primitive that cannot be represented by the component model.
    pub fn parse(text: &str) -> Result<Self, SmartsError> {
        crate::smarts_parse::parse(text)
    }

    /// Finds every distinct atom mapping in deterministic lexicographic order.
    #[must_use]
    pub fn find_matches(&self, component: &Component) -> Vec<SmartsMatch> {
        crate::smarts_match::find_matches(self, component)
    }

    /// Finds mappings in an annotated structure using its shared bond table.
    ///
    /// # Errors
    ///
    /// Returns [`SmartsDataError`] when connectivity or a chemical annotation
    /// required by this pattern is unavailable.
    pub fn find_structure_matches(
        &self,
        structure: &Structure,
    ) -> Result<Vec<SmartsMatch>, SmartsDataError> {
        crate::smarts_match::find_structure_matches(self, structure)
    }

    /// Reports whether at least one substructure mapping exists.
    #[must_use]
    pub fn matches(&self, component: &Component) -> bool {
        crate::smarts_match::has_match(self, component)
    }

    /// Number of atoms in the query graph.
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}

/// One SMARTS match, mapping query atom order to component atom indices.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SmartsMatch {
    /// Component atom index for every query atom.
    pub atom_indices: Box<[usize]>,
}

/// A syntax or unsupported-data error with a byte position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmartsError {
    /// Zero-based byte position in the input.
    pub position: usize,
    /// Human-readable failure reason.
    pub message: Box<str>,
}

impl SmartsError {
    pub(crate) fn new(position: usize, message: impl Into<Box<str>>) -> Self {
        Self {
            position,
            message: message.into(),
        }
    }
}

impl fmt::Display for SmartsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "SMARTS error at byte {}: {}",
            self.position, self.message
        )
    }
}

impl std::error::Error for SmartsError {}

/// Structure data required by a SMARTS primitive was not available.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmartsDataError {
    /// The structure has no resolved bond graph.
    Connectivity,
    /// Aromatic atom annotations are required.
    Aromaticity,
    /// CCD-backed formal charges are required.
    FormalCharge,
    /// Modelled stereochemical annotations are required.
    Stereochemistry,
}

impl fmt::Display for SmartsDataError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Connectivity => "resolved connectivity is unavailable",
            Self::Aromaticity => "aromatic atom annotations are unavailable",
            Self::FormalCharge => "formal-charge annotations are unavailable",
            Self::Stereochemistry => "stereochemical annotations are unavailable",
        })
    }
}

impl std::error::Error for SmartsDataError {}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AtomExpression {
    /// Comma-separated alternatives, each containing AND-connected tests.
    pub(crate) alternatives: Vec<Vec<SignedAtomTest>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SignedAtomTest {
    pub(crate) negated: bool,
    pub(crate) test: AtomTest,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum AtomTest {
    Any,
    Element(Element),
    Aromatic,
    Aliphatic,
    Degree(u8),
    Connectivity(u8),
    Hydrogens(u8),
    Charge(i8),
    RingCount(Option<u8>),
    RingSize(Option<u8>),
    Valence(u8),
    RingBonds(u8),
    Stereo(StereoConfiguration),
    Recursive(Box<SmartsPattern>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PatternBond {
    pub(crate) first: usize,
    pub(crate) second: usize,
    pub(crate) expression: BondExpression,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BondExpression {
    Default,
    Any,
    Order(BondOrder),
    Ring(bool),
    NotOrder(BondOrder),
}
