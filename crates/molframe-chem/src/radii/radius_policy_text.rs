//! The names a front end uses for the radius policy: one spelling, one parser.
//!
//! The command line and Python pass the same words, so the words are parsed
//! here, once, and printed back by `Display` for provenance records.

use crate::radius_policy::{HeteroPolicy, MetalFallback};
use core::fmt;
use core::str::FromStr;

/// A policy word that names nothing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolicyWordError {
    kind: &'static str,
    word: String,
    expected: &'static str,
}

impl fmt::Display for PolicyWordError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            out,
            "unknown {} `{}`; expected {}",
            self.kind, self.word, self.expected
        )
    }
}

impl core::error::Error for PolicyWordError {}

impl fmt::Display for HeteroPolicy {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if *self == Self::POLYMER {
            return out.write_str("none");
        }
        if *self == Self::ALL {
            return out.write_str("all");
        }
        if *self == Self::SOLUTE {
            return out.write_str("solute");
        }
        let classes = [
            (self.waters, "waters"),
            (self.ligands, "ligands"),
            (self.ions, "ions"),
        ];
        let named: Vec<&str> = classes
            .iter()
            .filter(|(on, _)| *on)
            .map(|(_, name)| *name)
            .collect();
        out.write_str(&named.join(","))
    }
}

impl FromStr for HeteroPolicy {
    type Err = PolicyWordError;

    /// Accepts `none` (the polymer alone), `solute` (ligands and ions), `all`,
    /// or a comma-separated union of `waters`, `ligands` and `ions`.
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let mut policy = Self::POLYMER;
        for word in text.split(',').map(str::trim) {
            match word {
                "none" | "polymer" | "" => {}
                "all" => policy = Self::ALL,
                "solute" => {
                    policy.ligands = true;
                    policy.ions = true;
                }
                "waters" | "water" => policy.waters = true,
                "ligands" | "ligand" => policy.ligands = true,
                "ions" | "ion" => policy.ions = true,
                other => {
                    return Err(PolicyWordError {
                        kind: "heteroatom class",
                        word: other.to_owned(),
                        expected: "none, solute, all, or a list of waters, ligands, ions",
                    });
                }
            }
        }
        Ok(policy)
    }
}

impl fmt::Display for MetalFallback {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.write_str(match self {
            Self::None => "none",
            Self::Vdw => "vdw",
            Self::Ionic => "ionic",
        })
    }
}

impl FromStr for MetalFallback {
    type Err = PolicyWordError;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        match text.trim() {
            "none" => Ok(Self::None),
            "vdw" => Ok(Self::Vdw),
            "ionic" => Ok(Self::Ionic),
            other => Err(PolicyWordError {
                kind: "metal fallback",
                word: other.to_owned(),
                expected: "none, vdw or ionic",
            }),
        }
    }
}

#[cfg(test)]
#[path = "radius_policy_text_tests.rs"]
mod tests;
