//! Density-map summary statistics and deterministic histograms.

use crate::DensityMap;
use crate::numeric::{f64_to_usize, usize_to_f64};

/// Summary over a selected set of finite density values.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MapStatistics {
    /// Number of included voxels.
    pub count: usize,
    /// Smallest included density.
    pub minimum: f32,
    /// Largest included density.
    pub maximum: f32,
    /// Arithmetic mean.
    pub mean: f64,
    /// Population standard deviation.
    pub sigma: f64,
}

/// Equally spaced density histogram.
#[derive(Clone, Debug, PartialEq)]
pub struct MapHistogram {
    /// Inclusive lower density bound.
    pub minimum: f32,
    /// Inclusive upper density bound; exact upper-bound values enter the last bin.
    pub maximum: f32,
    /// Count in each increasing density interval.
    pub counts: Vec<usize>,
}

/// Invalid statistics request or density field.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum MapStatisticsError {
    /// No values were selected.
    #[error("density-map statistics require at least one selected voxel")]
    Empty,
    /// The mask length does not match the map.
    #[error("density-map mask length does not match the grid")]
    MaskLength,
    /// A selected density is not finite.
    #[error("selected density-map voxel is non-finite")]
    NonFiniteDensity,
    /// Histogram bounds or bin count are invalid.
    #[error("density-map histogram range or bin count is invalid")]
    InvalidHistogram,
}

impl DensityMap {
    /// Calculates mean, population sigma, and extrema over all voxels.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty map or a non-finite density.
    pub fn statistics(&self) -> Result<MapStatistics, MapStatisticsError> {
        calculate(self.values.iter().copied())
    }

    /// Calculates statistics over voxels selected by a parallel mask.
    ///
    /// # Errors
    ///
    /// Returns an error for a mismatched mask, no selected values, or a
    /// selected non-finite density.
    pub fn masked_statistics(&self, mask: &[bool]) -> Result<MapStatistics, MapStatisticsError> {
        if mask.len() != self.values.len() {
            return Err(MapStatisticsError::MaskLength);
        }
        calculate(
            self.values
                .iter()
                .copied()
                .zip(mask)
                .filter_map(|(value, include)| include.then_some(value)),
        )
    }

    /// Builds an equally spaced histogram over an explicit density range.
    /// Values outside the range are omitted.
    ///
    /// # Errors
    ///
    /// Returns an error for zero bins, non-finite or reversed bounds, or a
    /// non-finite density value.
    pub fn histogram(
        &self,
        bins: usize,
        minimum: f32,
        maximum: f32,
    ) -> Result<MapHistogram, MapStatisticsError> {
        if bins == 0 || !minimum.is_finite() || !maximum.is_finite() || minimum >= maximum {
            return Err(MapStatisticsError::InvalidHistogram);
        }
        let mut counts = vec![0_usize; bins];
        let width = f64::from(maximum - minimum) / usize_to_f64(bins);
        for value in &self.values {
            if !value.is_finite() {
                return Err(MapStatisticsError::NonFiniteDensity);
            }
            if *value < minimum || *value > maximum {
                continue;
            }
            let scaled = f64_to_usize((f64::from(*value - minimum) / width).floor()).min(bins - 1);
            counts[scaled] += 1;
        }
        Ok(MapHistogram {
            minimum,
            maximum,
            counts,
        })
    }
}

fn calculate<I>(values: I) -> Result<MapStatistics, MapStatisticsError>
where
    I: IntoIterator<Item = f32>,
{
    let mut count = 0_usize;
    let mut minimum = f32::INFINITY;
    let mut maximum = f32::NEG_INFINITY;
    let mut mean = 0.0_f64;
    let mut sum_squares = 0.0_f64;
    for value in values {
        if !value.is_finite() {
            return Err(MapStatisticsError::NonFiniteDensity);
        }
        count += 1;
        minimum = minimum.min(value);
        maximum = maximum.max(value);
        let delta = f64::from(value) - mean;
        mean += delta / usize_to_f64(count);
        sum_squares += delta * (f64::from(value) - mean);
    }
    if count == 0 {
        return Err(MapStatisticsError::Empty);
    }
    Ok(MapStatistics {
        count,
        minimum,
        maximum,
        mean,
        sigma: (sum_squares / usize_to_f64(count)).sqrt(),
    })
}

#[cfg(test)]
#[path = "map_statistics_tests.rs"]
mod tests;
