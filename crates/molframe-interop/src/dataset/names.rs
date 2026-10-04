//! The words a split strategy is spelled with.

use super::{DatasetWarning, SplitStrategy};
use molframe_core::contract::{PolicyParseError, canonical_spelling};

impl SplitStrategy {
    /// The words a strategy can be spelled with, in declaration order.
    pub const NAMES: &'static [&'static str] = &[
        "sequence-identity",
        "structural-cluster",
        "temporal",
        "random",
    ];

    /// Reads a strategy from its name and the parameter it needs.
    ///
    /// `threshold` is required by `sequence-identity`, `seed` by `random`;
    /// `structural-cluster` and `temporal` take neither. A parameter the
    /// strategy does not use is refused rather than ignored, and a random split
    /// never invents a seed.
    ///
    /// # Errors
    ///
    /// Returns [`PolicyParseError`] for an unknown name, a missing parameter or
    /// a parameter the named strategy does not take.
    pub fn from_parts(
        name: &str,
        threshold: Option<f64>,
        seed: Option<u64>,
    ) -> Result<Self, PolicyParseError> {
        let refused = |expected: &str| PolicyParseError::new("strategy", name, expected);
        match (canonical_spelling(name).as_str(), threshold, seed) {
            ("sequence-identity", Some(threshold), None) => {
                Ok(Self::SequenceIdentity { threshold })
            }
            ("sequence-identity", _, _) => {
                Err(refused("sequence-identity with a threshold and no seed"))
            }
            ("structural-cluster", None, None) => Ok(Self::StructuralCluster),
            ("structural-cluster", _, _) => Err(refused(
                "structural-cluster, which takes no threshold or seed",
            )),
            ("temporal", None, None) => Ok(Self::Temporal),
            ("temporal", _, _) => Err(refused("temporal, which takes no threshold or seed")),
            ("random", None, Some(seed)) => Ok(Self::Random { seed }),
            ("random", _, _) => Err(refused("random with a seed and no threshold")),
            _ => Err(refused(&Self::NAMES.join(", "))),
        }
    }
}

impl std::fmt::Display for DatasetWarning {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RandomSplitMayLeak => formatter.write_str(
                "a random split makes no redundancy guarantee: related structures may land \
                 in different partitions",
            ),
        }
    }
}

#[cfg(test)]
#[path = "names_tests.rs"]
mod tests;
