//! Hausdorff and discrete Fréchet similarity between coordinate paths.

use pdbiox_geom::{rmsd, superpose};

/// Frame-to-frame distance used by path comparison.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PathFrameMetric {
    /// RMSD in the supplied Cartesian frame.
    CartesianRmsd,
    /// RMSD after optimal rigid superposition of every frame pair.
    FittedRmsd,
}

/// Two complementary geometric distances between paths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathSimilarity {
    /// Symmetric maximum-of-minimum frame distance.
    pub hausdorff_distance: f64,
    /// Order-sensitive discrete Fréchet distance.
    pub discrete_frechet_distance: f64,
}

/// Why path similarity could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum PathSimilarityError {
    /// Paths must contain finite, non-empty, equal-width coordinate frames.
    #[error("paths must contain finite, non-empty frames with a shared atom count")]
    InvalidPaths,
    /// The frame-distance matrix exceeds the caller's memory ceiling.
    #[error("path distance matrix exceeds the requested memory limit")]
    MemoryLimit,
    /// A fitted frame pair is geometrically degenerate.
    #[error("a requested rigid frame fit is degenerate")]
    DegenerateFit,
}

/// Compares two ordered coordinate paths under an explicit frame metric.
///
/// # Errors
///
/// Returns [`PathSimilarityError`] for invalid paths, allocation, or rigid fitting.
pub fn path_similarity(
    first: &[Vec<[f32; 3]>],
    second: &[Vec<[f32; 3]>],
    metric: PathFrameMetric,
    memory_limit_bytes: usize,
) -> Result<PathSimilarity, PathSimilarityError> {
    validate(first, second)?;
    let elements = first
        .len()
        .checked_mul(second.len())
        .ok_or(PathSimilarityError::MemoryLimit)?;
    let required = elements
        .checked_mul(size_of::<f64>())
        .ok_or(PathSimilarityError::MemoryLimit)?;
    if required > memory_limit_bytes {
        return Err(PathSimilarityError::MemoryLimit);
    }
    let mut distances = vec![0.0; elements];
    for (row, first_frame) in first.iter().enumerate() {
        for (column, second_frame) in second.iter().enumerate() {
            distances[row * second.len() + column] =
                frame_distance(first_frame, second_frame, metric)?;
        }
    }
    Ok(PathSimilarity {
        hausdorff_distance: hausdorff(&distances, first.len(), second.len()),
        discrete_frechet_distance: discrete_frechet(&distances, first.len(), second.len()),
    })
}

fn validate(first: &[Vec<[f32; 3]>], second: &[Vec<[f32; 3]>]) -> Result<(), PathSimilarityError> {
    let Some(atom_count) = first.first().map(Vec::len) else {
        return Err(PathSimilarityError::InvalidPaths);
    };
    if atom_count == 0
        || second.is_empty()
        || first.iter().chain(second).any(|frame| {
            frame.len() != atom_count || frame.iter().flatten().any(|value| !value.is_finite())
        })
    {
        return Err(PathSimilarityError::InvalidPaths);
    }
    Ok(())
}

fn frame_distance(
    first: &[[f32; 3]],
    second: &[[f32; 3]],
    metric: PathFrameMetric,
) -> Result<f64, PathSimilarityError> {
    match metric {
        PathFrameMetric::CartesianRmsd => {
            rmsd(first, second).map_err(|_| PathSimilarityError::InvalidPaths)
        }
        PathFrameMetric::FittedRmsd => superpose(first, second)
            .map(|fit| fit.rmsd)
            .map_err(|_| PathSimilarityError::DegenerateFit),
    }
}

fn hausdorff(distances: &[f64], rows: usize, columns: usize) -> f64 {
    let first_to_second = (0..rows)
        .map(|row| {
            (0..columns)
                .map(|column| distances[row * columns + column])
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    let second_to_first = (0..columns)
        .map(|column| {
            (0..rows)
                .map(|row| distances[row * columns + column])
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    first_to_second.max(second_to_first)
}

fn discrete_frechet(distances: &[f64], rows: usize, columns: usize) -> f64 {
    let mut previous: Vec<f64> = vec![0.0; columns];
    let mut current: Vec<f64> = vec![0.0; columns];
    for row in 0..rows {
        for column in 0..columns {
            let distance = distances[row * columns + column];
            current[column] = match (row, column) {
                (0, 0) => distance,
                (0, _) => current[column - 1].max(distance),
                (_, 0) => previous[0].max(distance),
                _ => previous[column]
                    .min(previous[column - 1])
                    .min(current[column - 1])
                    .max(distance),
            };
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[columns - 1]
}

#[cfg(test)]
#[path = "path_similarity_tests.rs"]
mod tests;
