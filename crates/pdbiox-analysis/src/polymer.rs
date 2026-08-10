//! Polymer contour, end-to-end, and tangent-correlation statistics.

use crate::numeric::usize_to_f64;

/// Geometric statistics for one explicitly ordered polymer path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolymerStatistics {
    /// Sum of consecutive segment lengths.
    pub contour_length: f64,
    /// Distance between the first and last path sites.
    pub end_to_end_distance: f64,
    /// Mean consecutive segment length.
    pub mean_segment_length: f64,
    /// Mean dot product of adjacent unit tangents.
    pub adjacent_tangent_correlation: Option<f64>,
    /// Discrete worm-like-chain estimate `-mean_segment / ln(correlation)`.
    pub persistence_length: Option<f64>,
}

/// Why polymer statistics could not be computed.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum PolymerError {
    /// A polymer path needs at least two sites.
    #[error("a polymer path requires at least two ordered sites")]
    TooShort,
    /// Coordinates must be finite.
    #[error("polymer coordinates must be finite")]
    NonFinite,
    /// Consecutive sites must not coincide.
    #[error("consecutive polymer sites must not coincide")]
    CoincidentSites,
}

/// Computes statistics for an explicitly ordered polymer path.
///
/// The path may represent backbone atoms, monomer centres, or another caller
/// definition. No residue names or canonical atom names are inferred. The
/// persistence estimate is present only when at least two segments yield a
/// positive correlation strictly below one; other geometries are reported as
/// `None` instead of receiving a fabricated finite value.
///
/// # Errors
///
/// Returns [`PolymerError`] for short, non-finite, or degenerate paths.
pub fn polymer_statistics(path: &[[f32; 3]]) -> Result<PolymerStatistics, PolymerError> {
    if path.len() < 2 {
        return Err(PolymerError::TooShort);
    }
    if path.iter().flatten().any(|value| !value.is_finite()) {
        return Err(PolymerError::NonFinite);
    }
    let mut lengths = Vec::with_capacity(path.len() - 1);
    let mut tangents = Vec::with_capacity(path.len() - 1);
    for pair in path.windows(2) {
        let displacement = subtract(pair[1], pair[0]);
        let length = norm(displacement);
        if length == 0.0 {
            return Err(PolymerError::CoincidentSites);
        }
        lengths.push(length);
        tangents.push(displacement.map(|value| value / length));
    }
    let contour_length: f64 = lengths.iter().sum();
    let mean_segment_length = contour_length / usize_to_f64(lengths.len());
    let adjacent_tangent_correlation = if tangents.len() < 2 {
        None
    } else {
        Some(
            tangents
                .windows(2)
                .map(|pair| dot(pair[0], pair[1]))
                .sum::<f64>()
                / usize_to_f64(tangents.len() - 1),
        )
    };
    let persistence_length = adjacent_tangent_correlation
        .filter(|correlation| *correlation > 0.0 && *correlation < 1.0)
        .map(|correlation| -mean_segment_length / correlation.ln());
    Ok(PolymerStatistics {
        contour_length,
        end_to_end_distance: norm(subtract(path[path.len() - 1], path[0])),
        mean_segment_length,
        adjacent_tangent_correlation,
        persistence_length,
    })
}

fn subtract(left: [f32; 3], right: [f32; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| f64::from(left[axis] - right[axis]))
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left.into_iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

#[cfg(test)]
#[path = "polymer_tests.rs"]
mod tests;
