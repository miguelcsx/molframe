//! The words a caller uses to name an ensemble method, read once, here.

use crate::{Linkage, PathFrameMetric, RemainderPolicy};
use molframe_core::{Code, diagnostic_from};
use std::str::FromStr;

/// A word that names none of the choices of its decision.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{word:?} is not a {decision}; the choices are {choices}")]
pub struct UnknownWord {
    decision: &'static str,
    word: Box<str>,
    choices: &'static str,
}

diagnostic_from!(UnknownWord, |_error| Code::E5101);

impl FromStr for Linkage {
    type Err = UnknownWord;

    fn from_str(word: &str) -> Result<Self, Self::Err> {
        match word {
            "single" => Ok(Self::Single),
            "complete" => Ok(Self::Complete),
            "average" => Ok(Self::Average),
            other => Err(UnknownWord {
                decision: "linkage",
                word: other.into(),
                choices: "single, complete and average",
            }),
        }
    }
}

impl FromStr for RemainderPolicy {
    type Err = UnknownWord;

    fn from_str(word: &str) -> Result<Self, Self::Err> {
        match word {
            "include" => Ok(Self::Include),
            "reject" => Ok(Self::Reject),
            other => Err(UnknownWord {
                decision: "remainder policy",
                word: other.into(),
                choices: "include and reject",
            }),
        }
    }
}

impl FromStr for PathFrameMetric {
    type Err = UnknownWord;

    fn from_str(word: &str) -> Result<Self, Self::Err> {
        match word {
            "cartesian_rmsd" | "cartesian-rmsd" => Ok(Self::CartesianRmsd),
            "fitted_rmsd" | "fitted-rmsd" => Ok(Self::FittedRmsd),
            other => Err(UnknownWord {
                decision: "path frame metric",
                word: other.into(),
                choices: "cartesian_rmsd and fitted_rmsd",
            }),
        }
    }
}

#[cfg(test)]
#[path = "words_tests.rs"]
mod tests;
