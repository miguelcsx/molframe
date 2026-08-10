//! Circular angles and product-torus distances.
//!
//! A torsion is not a scalar on a line: values separated by a whole turn are
//! identical. These types keep that topology explicit and make every distance
//! use the shortest wrapped displacement. Construction and every summary are
//! deterministic and reject non-finite inputs.

use crate::numeric::exact_count;
use std::f64::consts::{PI, TAU};

/// Failure to construct or compare periodic angular data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodicError {
    /// An angle, tolerance or weight was not finite.
    NonFinite,
    /// Two points or a point and its metric have different dimensions.
    DimensionMismatch,
    /// A metric weight was zero or negative.
    InvalidWeight,
    /// No observation was supplied.
    Empty,
    /// The observation count could not be represented by the numeric kernel.
    TooManyObservations,
}

/// One finite angle canonically stored in `(-π, π]`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct PeriodicAngle(f64);

impl PeriodicAngle {
    /// Builds an angle from radians.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::NonFinite`] for NaN or infinity.
    pub fn from_radians(radians: f64) -> Result<Self, PeriodicError> {
        if !radians.is_finite() {
            return Err(PeriodicError::NonFinite);
        }
        Ok(Self(wrap(radians)))
    }

    /// Builds an angle from degrees.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::NonFinite`] for NaN or infinity.
    pub fn from_degrees(degrees: f64) -> Result<Self, PeriodicError> {
        Self::from_radians(degrees.to_radians())
    }

    /// The canonical value in radians.
    #[must_use]
    pub const fn radians(self) -> f64 {
        self.0
    }

    /// The canonical value in degrees.
    #[must_use]
    pub fn degrees(self) -> f64 {
        self.0.to_degrees()
    }

    /// The shortest signed displacement carrying `self` to `other`.
    #[must_use]
    pub fn signed_delta(self, other: Self) -> f64 {
        wrap(other.0 - self.0)
    }

    /// The unsigned geodesic distance on the circle.
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        self.signed_delta(other).abs()
    }
}

/// Circular location and dispersion for a sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularSummary {
    /// Mean direction, absent when the sample has no preferred direction.
    pub mean: Option<PeriodicAngle>,
    /// Length of the mean resultant vector in `[0, 1]`.
    pub resultant: f64,
    /// Circular variance, `1 - resultant`.
    pub variance: f64,
}

/// Computes the circular mean and dispersion.
///
/// `indeterminate_tolerance` is the largest resultant length treated as having
/// no defensible mean direction.
///
/// # Errors
///
/// Returns an error for an empty sample or an invalid tolerance.
pub fn circular_summary(
    angles: &[PeriodicAngle],
    indeterminate_tolerance: f64,
) -> Result<CircularSummary, PeriodicError> {
    if angles.is_empty() {
        return Err(PeriodicError::Empty);
    }
    if !indeterminate_tolerance.is_finite() || indeterminate_tolerance < 0.0 {
        return Err(PeriodicError::NonFinite);
    }
    let (sin_sum, cos_sum) = angles.iter().fold((0.0, 0.0), |(sin, cos), angle| {
        (sin + angle.0.sin(), cos + angle.0.cos())
    });
    let count = exact_count(angles.len()).ok_or(PeriodicError::TooManyObservations)?;
    let resultant = sin_sum.hypot(cos_sum) / count;
    let mean = if resultant <= indeterminate_tolerance {
        None
    } else {
        Some(PeriodicAngle(wrap(sin_sum.atan2(cos_sum))))
    };
    Ok(CircularSummary {
        mean,
        resultant,
        variance: 1.0 - resultant,
    })
}

/// Weighted product metric over a fixed number of circles.
#[derive(Clone, Debug, PartialEq)]
pub struct TorusMetric {
    weights: Box<[f64]>,
}

impl TorusMetric {
    /// Builds a metric from strictly positive finite weights.
    ///
    /// # Errors
    ///
    /// Refuses empty, non-finite or non-positive weights.
    pub fn new(weights: impl Into<Box<[f64]>>) -> Result<Self, PeriodicError> {
        let weights = weights.into();
        if weights.is_empty() {
            return Err(PeriodicError::Empty);
        }
        if weights.iter().any(|weight| !weight.is_finite()) {
            return Err(PeriodicError::NonFinite);
        }
        if weights.iter().any(|weight| *weight <= 0.0) {
            return Err(PeriodicError::InvalidWeight);
        }
        Ok(Self { weights })
    }

    /// Builds the unweighted metric in `dimension` dimensions.
    ///
    /// # Errors
    ///
    /// Refuses a zero-dimensional torus.
    pub fn uniform(dimension: usize) -> Result<Self, PeriodicError> {
        Self::new(vec![1.0; dimension].into_boxed_slice())
    }

    /// The number of angular coordinates expected by this metric.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.weights.len()
    }

    /// Computes the weighted geodesic distance on the product torus.
    ///
    /// # Errors
    ///
    /// Refuses dimensions that do not match the metric.
    pub fn distance(
        &self,
        left: &[PeriodicAngle],
        right: &[PeriodicAngle],
    ) -> Result<f64, PeriodicError> {
        if left.len() != self.weights.len() || right.len() != self.weights.len() {
            return Err(PeriodicError::DimensionMismatch);
        }
        let squared = left
            .iter()
            .zip(right)
            .zip(&self.weights)
            .map(|((left, right), weight)| weight * left.signed_delta(*right).powi(2))
            .sum::<f64>();
        Ok(squared.sqrt())
    }
}

/// Component-wise Fréchet summaries on a product torus.
///
/// # Errors
///
/// Refuses an empty sample, ragged points or an invalid tolerance.
pub fn torus_summary(
    points: &[Box<[PeriodicAngle]>],
    indeterminate_tolerance: f64,
) -> Result<Box<[CircularSummary]>, PeriodicError> {
    let Some(first) = points.first() else {
        return Err(PeriodicError::Empty);
    };
    if first.is_empty() || points.iter().any(|point| point.len() != first.len()) {
        return Err(PeriodicError::DimensionMismatch);
    }
    (0..first.len())
        .map(|dimension| {
            let values: Vec<_> = points.iter().map(|point| point[dimension]).collect();
            circular_summary(&values, indeterminate_tolerance)
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Vec::into_boxed_slice)
}

fn wrap(radians: f64) -> f64 {
    let wrapped = (radians + PI).rem_euclid(TAU) - PI;
    if wrapped <= -PI { PI } else { wrapped }
}

#[cfg(test)]
#[path = "periodic_angle_tests.rs"]
mod tests;
