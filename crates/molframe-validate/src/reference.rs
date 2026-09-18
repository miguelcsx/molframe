//! Versioned empirical reference distributions for validation observations.

use crate::{BondDeviation, RamachandranRecord};
use std::collections::BTreeMap;

/// A validated one- or two-dimensional empirical histogram.
#[derive(Clone, Debug, PartialEq)]
pub enum ReferenceDistribution {
    /// Scalar bins, such as signed bond-length deviations.
    Histogram {
        /// Stable metric name inside a reference set.
        name: Box<str>,
        /// Strictly increasing bin boundaries.
        edges: Vec<f64>,
        /// Non-negative empirical weight per interval.
        weights: Vec<f64>,
        /// Cached total mass, validated as finite and positive.
        total: f64,
    },
    /// Rectangular grid, such as a φ/ψ distribution.
    Grid {
        /// Stable metric name inside a reference set.
        name: Box<str>,
        /// Strictly increasing first-axis boundaries.
        x_edges: Vec<f64>,
        /// Strictly increasing second-axis boundaries.
        y_edges: Vec<f64>,
        /// Row-major empirical weights: X interval first, then Y.
        weights: Vec<f64>,
        /// Cached total mass, validated as finite and positive.
        total: f64,
    },
}

impl ReferenceDistribution {
    /// Builds a validated scalar histogram.
    ///
    /// # Errors
    ///
    /// Refuses an empty name, invalid edges, mismatched weights or zero mass.
    pub fn histogram(
        name: impl Into<Box<str>>,
        edges: Vec<f64>,
        weights: Vec<f64>,
    ) -> Result<Self, ReferenceError> {
        let name = name.into();
        validate_name(&name)?;
        validate_axis(&edges)?;
        if weights.len() + 1 != edges.len() {
            return Err(ReferenceError::Shape);
        }
        let total = validate_weights(&weights)?;
        Ok(Self::Histogram {
            name,
            edges,
            weights,
            total,
        })
    }

    /// Builds a validated two-dimensional empirical grid.
    ///
    /// # Errors
    ///
    /// Refuses an empty name, invalid axes, mismatched weights or zero mass.
    pub fn grid(
        name: impl Into<Box<str>>,
        x_edges: Vec<f64>,
        y_edges: Vec<f64>,
        weights: Vec<f64>,
    ) -> Result<Self, ReferenceError> {
        let name = name.into();
        validate_name(&name)?;
        validate_axis(&x_edges)?;
        validate_axis(&y_edges)?;
        let expected = (x_edges.len() - 1)
            .checked_mul(y_edges.len() - 1)
            .ok_or(ReferenceError::Shape)?;
        if weights.len() != expected {
            return Err(ReferenceError::Shape);
        }
        let total = validate_weights(&weights)?;
        Ok(Self::Grid {
            name,
            x_edges,
            y_edges,
            weights,
            total,
        })
    }

    /// Stable name used to select this distribution.
    #[must_use]
    pub fn name(&self) -> &str {
        match self {
            Self::Histogram { name, .. } | Self::Grid { name, .. } => name,
        }
    }
}

/// Immutable collection whose identity includes the exact reference version.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceLibrary {
    id: Box<str>,
    version: Box<str>,
    distributions: BTreeMap<Box<str>, ReferenceDistribution>,
}

impl ReferenceLibrary {
    /// Creates a named, versioned set and rejects duplicate metric names.
    ///
    /// # Errors
    ///
    /// Refuses empty identity fields, no distributions or duplicate names.
    pub fn new(
        id: impl Into<Box<str>>,
        version: impl Into<Box<str>>,
        distributions: impl IntoIterator<Item = ReferenceDistribution>,
    ) -> Result<Self, ReferenceError> {
        let id = id.into();
        let version = version.into();
        validate_name(&id)?;
        validate_name(&version)?;
        let mut indexed = BTreeMap::new();
        for distribution in distributions {
            let name: Box<str> = distribution.name().into();
            if indexed.insert(name, distribution).is_some() {
                return Err(ReferenceError::Duplicate);
            }
        }
        if indexed.is_empty() {
            return Err(ReferenceError::Empty);
        }
        Ok(Self {
            id,
            version,
            distributions: indexed,
        })
    }

    /// Collection identifier independent of its release.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Exact reference-data release.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Verifies that a name resolves to a two-dimensional grid.
    ///
    /// This lets an analysis validate its configuration once, before walking
    /// observations, without exposing the library's storage representation.
    ///
    /// # Errors
    ///
    /// Returns missing or dimension errors for an absent name or histogram.
    pub fn validate_grid(&self, name: &str) -> Result<(), ReferenceError> {
        match self.distributions.get(name) {
            Some(ReferenceDistribution::Grid { .. }) => Ok(()),
            Some(ReferenceDistribution::Histogram { .. }) => Err(ReferenceError::Dimension),
            None => Err(ReferenceError::Missing),
        }
    }

    /// Assesses one scalar against a named histogram.
    ///
    /// # Errors
    ///
    /// Returns missing, dimension or non-finite observation errors.
    pub fn assess_scalar(
        &self,
        name: &str,
        value: f64,
    ) -> Result<ReferenceAssessment, ReferenceError> {
        if !value.is_finite() {
            return Err(ReferenceError::Observation);
        }
        let Some(ReferenceDistribution::Histogram {
            edges,
            weights,
            total,
            ..
        }) = self.distributions.get(name)
        else {
            return self.missing_or_dimension(name);
        };
        let Some(bin) = bin(edges, value) else {
            return Ok(self.assessment(name, 0.0, Some(if value < edges[0] { 0.0 } else { 1.0 })));
        };
        let preceding: f64 = weights[..bin].iter().sum();
        let fraction = (value - edges[bin]) / (edges[bin + 1] - edges[bin]);
        Ok(self.assessment(
            name,
            weights[bin] / *total,
            Some((preceding + fraction * weights[bin]) / *total),
        ))
    }

    /// Assesses a point against a named two-dimensional grid.
    ///
    /// # Errors
    ///
    /// Returns missing, dimension or non-finite observation errors.
    pub fn assess_pair(
        &self,
        name: &str,
        x: f64,
        y: f64,
    ) -> Result<ReferenceAssessment, ReferenceError> {
        if !x.is_finite() || !y.is_finite() {
            return Err(ReferenceError::Observation);
        }
        let Some(ReferenceDistribution::Grid {
            x_edges,
            y_edges,
            weights,
            total,
            ..
        }) = self.distributions.get(name)
        else {
            return self.missing_or_dimension(name);
        };
        let probability = match (bin(x_edges, x), bin(y_edges, y)) {
            (Some(x_bin), Some(y_bin)) => weights[x_bin * (y_edges.len() - 1) + y_bin] / *total,
            _ => 0.0,
        };
        Ok(self.assessment(name, probability, None))
    }

    fn assessment(
        &self,
        distribution: &str,
        probability: f64,
        percentile: Option<f64>,
    ) -> ReferenceAssessment {
        ReferenceAssessment {
            set: self.id.clone(),
            version: self.version.clone(),
            distribution: distribution.into(),
            probability,
            percentile,
        }
    }

    fn missing_or_dimension<T>(&self, name: &str) -> Result<T, ReferenceError> {
        if self.distributions.contains_key(name) {
            Err(ReferenceError::Dimension)
        } else {
            Err(ReferenceError::Missing)
        }
    }
}

/// Reference identity and normalized empirical support for one observation.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceAssessment {
    /// Reference collection identifier.
    pub set: Box<str>,
    /// Exact collection version.
    pub version: Box<str>,
    /// Distribution name.
    pub distribution: Box<str>,
    /// Normalized mass of the containing bin, or zero outside the domain.
    pub probability: f64,
    /// Interpolated cumulative percentile for scalar histograms.
    pub percentile: Option<f64>,
}

/// Invalid reference data or observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ReferenceError {
    /// Identity or distribution name is empty.
    #[error("reference identity fields must not be empty")]
    Empty,
    /// Axis boundaries are not finite and strictly increasing.
    #[error("reference axis must be finite and strictly increasing")]
    Axis,
    /// Weight array shape does not match the declared axes.
    #[error("reference weight shape does not match its axes")]
    Shape,
    /// Weights are invalid or have no positive mass.
    #[error("reference weights must be finite, non-negative and have positive mass")]
    Weight,
    /// Two distributions use the same name.
    #[error("reference distribution names must be unique")]
    Duplicate,
    /// Requested distribution does not exist.
    #[error("reference distribution is absent")]
    Missing,
    /// Requested assessment has the wrong dimensionality.
    #[error("reference distribution dimensionality does not match the observation")]
    Dimension,
    /// Observation is not finite.
    #[error("reference observation must be finite")]
    Observation,
}

/// Assesses a measured bond deviation against a scalar reference histogram.
///
/// # Errors
///
/// Returns the errors documented by [`ReferenceLibrary::assess_scalar`].
pub fn assess_bond_deviation(
    deviation: &BondDeviation,
    references: &ReferenceLibrary,
    distribution: &str,
) -> Result<ReferenceAssessment, ReferenceError> {
    references.assess_scalar(distribution, f64::from(deviation.deviation))
}

/// Assesses a φ/ψ observation against a two-dimensional reference grid.
///
/// # Errors
///
/// Returns the errors documented by [`ReferenceLibrary::assess_pair`].
pub fn assess_ramachandran(
    record: &RamachandranRecord,
    references: &ReferenceLibrary,
    distribution: &str,
) -> Result<ReferenceAssessment, ReferenceError> {
    references.assess_pair(distribution, record.phi, record.psi)
}

fn validate_name(name: &str) -> Result<(), ReferenceError> {
    if name.trim().is_empty() {
        Err(ReferenceError::Empty)
    } else {
        Ok(())
    }
}

fn validate_axis(edges: &[f64]) -> Result<(), ReferenceError> {
    if edges.len() < 2
        || edges.iter().any(|edge| !edge.is_finite())
        || edges.windows(2).any(|pair| pair[0] >= pair[1])
    {
        Err(ReferenceError::Axis)
    } else {
        Ok(())
    }
}

fn validate_weights(weights: &[f64]) -> Result<f64, ReferenceError> {
    let total = weights.iter().sum::<f64>();
    if weights.is_empty()
        || weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight < 0.0)
        || !total.is_finite()
        || total <= 0.0
    {
        Err(ReferenceError::Weight)
    } else {
        Ok(total)
    }
}

fn bin(edges: &[f64], value: f64) -> Option<usize> {
    if value < edges[0] || value > edges[edges.len() - 1] {
        return None;
    }
    if value >= edges[edges.len() - 1] {
        return Some(edges.len() - 2);
    }
    edges.partition_point(|edge| *edge <= value).checked_sub(1)
}

#[cfg(test)]
#[path = "reference_tests.rs"]
mod tests;
