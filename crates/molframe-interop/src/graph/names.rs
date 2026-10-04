//! The words a graph request is spelled with.
//!
//! Names are kebab-case and an underscore is read as a hyphen, as in every
//! other vocabulary (`molframe_core::contract::PolicyParseError`).

use super::{EdgeDirection, EdgeFeature, EdgeKind, NodeFeature, NodeLevel};
use molframe_core::contract::{PolicyParseError, canonical_spelling};
use std::fmt;
use std::str::FromStr;

macro_rules! named {
    ($type:ty, $field:literal, [$(($variant:path, $name:literal)),+ $(,)?]) => {
        impl $type {
            /// The words this choice can be spelled with, in declaration order.
            pub const NAMES: &'static [&'static str] = &[$($name),+];

            /// The canonical spelling of this choice.
            #[must_use]
            pub const fn name(self) -> &'static str {
                match self {
                    $($variant => $name,)+
                }
            }
        }

        impl fmt::Display for $type {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(self.name())
            }
        }

        impl FromStr for $type {
            type Err = PolicyParseError;

            fn from_str(value: &str) -> Result<Self, Self::Err> {
                match canonical_spelling(value).as_str() {
                    $($name => Ok($variant),)+
                    _ => Err(PolicyParseError::new($field, value, &Self::NAMES.join(", "))),
                }
            }
        }
    };
}

named!(
    NodeLevel,
    "nodes",
    [
        (NodeLevel::Atoms, "atoms"),
        (NodeLevel::Residues, "residues")
    ]
);

named!(
    EdgeDirection,
    "direction",
    [
        (EdgeDirection::Undirected, "undirected"),
        (EdgeDirection::Symmetric, "symmetric"),
        (EdgeDirection::Directed, "directed"),
    ]
);

named!(
    NodeFeature,
    "node_feature",
    [
        (NodeFeature::Element, "element"),
        (NodeFeature::FormalCharge, "formal-charge"),
        (NodeFeature::PartialCharge, "partial-charge"),
        (NodeFeature::BFactor, "b-factor"),
        (NodeFeature::Occupancy, "occupancy"),
        (NodeFeature::AtomCount, "atom-count"),
        (NodeFeature::PositionX, "position-x"),
        (NodeFeature::PositionY, "position-y"),
        (NodeFeature::PositionZ, "position-z"),
    ]
);

named!(
    EdgeFeature,
    "edge_feature",
    [
        (EdgeFeature::Distance, "distance"),
        (EdgeFeature::BondOrder, "bond-order"),
    ]
);

impl EdgeKind {
    /// The words an edge algorithm can be spelled with.
    pub const NAMES: &'static [&'static str] = &["bonds", "contacts", "radius", "k-nearest"];

    /// Reads an algorithm from its name and the parameter it needs.
    ///
    /// `cutoff` (ångström) is required by `contacts` and `radius`, `neighbors`
    /// by `k-nearest`, and `bonds` takes neither. A parameter the algorithm
    /// does not use is refused rather than ignored.
    ///
    /// # Errors
    ///
    /// Returns [`PolicyParseError`] for an unknown name, a missing parameter or
    /// a parameter the named algorithm does not take.
    pub fn from_parts(
        name: &str,
        cutoff: Option<f32>,
        neighbors: Option<usize>,
    ) -> Result<Self, PolicyParseError> {
        let refused = |expected: &str| PolicyParseError::new("edges", name, expected);
        match (canonical_spelling(name).as_str(), cutoff, neighbors) {
            ("bonds", None, None) => Ok(Self::Bonds),
            ("bonds", _, _) => Err(refused("bonds, which takes no cutoff or neighbors")),
            ("contacts", Some(cutoff), None) => Ok(Self::Contacts { cutoff }),
            ("contacts", _, _) => Err(refused("contacts with a cutoff and no neighbors")),
            ("radius", Some(cutoff), None) => Ok(Self::Radius { cutoff }),
            ("radius", _, _) => Err(refused("radius with a cutoff and no neighbors")),
            ("k-nearest", None, Some(neighbors)) => Ok(Self::KNearest { neighbors }),
            ("k-nearest", _, _) => Err(refused("k-nearest with neighbors and no cutoff")),
            _ => Err(refused(&Self::NAMES.join(", "))),
        }
    }
}

#[cfg(test)]
#[path = "names_tests.rs"]
mod tests;
