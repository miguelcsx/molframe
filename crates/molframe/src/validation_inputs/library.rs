//! Versioned empirical reference distributions.

use super::ValidationInputError;
use crate::document::read_document;
use molframe_validate::{ReferenceDistribution, ReferenceLibrary};
use serde::Deserialize;
use std::path::Path;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    library: Identity,
    distribution: Vec<Distribution>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Identity {
    id: Box<str>,
    version: Box<str>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Distribution {
    name: Box<str>,
    edges: Option<Vec<f64>>,
    x_edges: Option<Vec<f64>>,
    y_edges: Option<Vec<f64>>,
    weights: Vec<f64>,
}

impl Distribution {
    /// A distribution is a scalar histogram (`edges`) or a grid (`x_edges` and
    /// `y_edges`); naming both forms, or half of the grid, is ambiguous.
    fn build(self) -> Result<ReferenceDistribution, ValidationInputError> {
        match (self.edges, self.x_edges, self.y_edges) {
            (Some(edges), None, None) => Ok(ReferenceDistribution::histogram(
                self.name,
                edges,
                self.weights,
            )?),
            (None, Some(x_edges), Some(y_edges)) => Ok(ReferenceDistribution::grid(
                self.name,
                x_edges,
                y_edges,
                self.weights,
            )?),
            _ => Err(ValidationInputError::invalid(
                format!("distribution {}", self.name),
                "give either `edges` for a histogram or both `x_edges` and `y_edges` for a grid",
            )),
        }
    }
}

/// Reads a versioned reference library from JSON or TOML.
///
/// A `[library]` table names the set (`id`, `version`), and each
/// `[[distribution]]` is either a scalar histogram with `name`, `edges`
/// (strictly increasing boundaries) and `weights` (one non-negative weight per
/// interval), or a two-dimensional grid with `name`, `x_edges`, `y_edges` and
/// `weights` in row-major order (first axis outermost). Keys equal the
/// parameters of [`ReferenceDistribution::histogram`] and
/// [`ReferenceDistribution::grid`]. Unknown keys are errors.
///
/// ```toml
/// [library]
/// id = "bond-deviation"
/// version = "2026-10"
///
/// [[distribution]]
/// name = "c-n"
/// edges = [-0.1, 0.0, 0.1]
/// weights = [1.0, 3.0]
/// ```
///
/// # Errors
///
/// Returns [`ValidationInputError`] for I/O, syntax or schema errors, an
/// ambiguous distribution, or anything [`ReferenceLibrary::new`] refuses.
pub fn read_reference_library(
    path: impl AsRef<Path>,
) -> Result<ReferenceLibrary, ValidationInputError> {
    let document: Document = read_document(path.as_ref())?;
    let distributions = document
        .distribution
        .into_iter()
        .map(Distribution::build)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ReferenceLibrary::new(
        document.library.id,
        document.library.version,
        distributions,
    )?)
}

#[cfg(test)]
#[path = "library_tests.rs"]
mod tests;
