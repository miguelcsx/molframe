//! Hausdorff and discrete Fréchet similarity between coordinate paths.
//!
//! Frame distances are consumed in row-major order instead of materialising a
//! pairwise matrix. Comparing paths with `F1` and `F2` frames of `A` atoms takes
//! `O(F1 * F2 * A)` time and `O(F2)` working memory.

use crate::frame_view::{FrameSource, FrameView};
use molframe_geom::{rmsd, superpose};

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
    path_similarity_source(first, second, metric, memory_limit_bytes)
}

/// Compares two borrowed contiguous coordinate paths without repacking frames.
///
/// # Errors
///
/// Returns [`PathSimilarityError`] for invalid paths, allocation, or rigid fitting.
pub fn path_similarity_view(
    first: FrameView<'_>,
    second: FrameView<'_>,
    metric: PathFrameMetric,
    memory_limit_bytes: usize,
) -> Result<PathSimilarity, PathSimilarityError> {
    path_similarity_source(&first, &second, metric, memory_limit_bytes)
}

fn path_similarity_source<F: FrameSource + ?Sized, S: FrameSource + ?Sized>(
    first: &F,
    second: &S,
    metric: PathFrameMetric,
    memory_limit_bytes: usize,
) -> Result<PathSimilarity, PathSimilarityError> {
    validate(first, second)?;
    let columns = second.frame_count();
    let workspace_elements = columns
        .checked_mul(2)
        .ok_or(PathSimilarityError::MemoryLimit)?;
    let required = workspace_elements
        .checked_mul(size_of::<f64>())
        .ok_or(PathSimilarityError::MemoryLimit)?;
    if required > memory_limit_bytes {
        return Err(PathSimilarityError::MemoryLimit);
    }
    let mut workspace = Vec::<f64>::new();
    workspace
        .try_reserve_exact(workspace_elements)
        .map_err(|_| PathSimilarityError::MemoryLimit)?;
    workspace.resize(workspace_elements, 0.0);
    let (frechet_row, column_minima) = workspace.split_at_mut(columns);
    column_minima.fill(f64::INFINITY);
    let mut first_to_second = 0.0_f64;
    for row in 0..first.frame_count() {
        let first_frame = first.frame(row).ok_or(PathSimilarityError::InvalidPaths)?;
        let mut row_minimum = f64::INFINITY;
        let mut previous_diagonal = 0.0;
        for column in 0..columns {
            let second_frame = second
                .frame(column)
                .ok_or(PathSimilarityError::InvalidPaths)?;
            let distance = frame_distance(first_frame, second_frame, metric)?;
            row_minimum = row_minimum.min(distance);
            column_minima[column] = column_minima[column].min(distance);

            let previous_row_value = frechet_row[column];
            frechet_row[column] = match (row, column) {
                (0, 0) => distance,
                (0, _) => frechet_row[column - 1].max(distance),
                (_, 0) => previous_row_value.max(distance),
                _ => previous_row_value
                    .min(previous_diagonal)
                    .min(frechet_row[column - 1])
                    .max(distance),
            };
            previous_diagonal = previous_row_value;
        }
        first_to_second = first_to_second.max(row_minimum);
    }
    let second_to_first = column_minima.iter().copied().fold(0.0, f64::max);
    Ok(PathSimilarity {
        hausdorff_distance: first_to_second.max(second_to_first),
        discrete_frechet_distance: frechet_row[columns - 1],
    })
}

fn validate<F: FrameSource + ?Sized, S: FrameSource + ?Sized>(
    first: &F,
    second: &S,
) -> Result<(), PathSimilarityError> {
    if first.frame_count() == 0 || second.frame_count() == 0 {
        return Err(PathSimilarityError::InvalidPaths);
    }
    let atom_count = first.atom_count();
    if atom_count == 0 || second.atom_count() != atom_count {
        return Err(PathSimilarityError::InvalidPaths);
    }
    for index in 0..first.frame_count() {
        let frame = first
            .frame(index)
            .ok_or(PathSimilarityError::InvalidPaths)?;
        if frame.len() != atom_count || frame.iter().flatten().any(|value| !value.is_finite()) {
            return Err(PathSimilarityError::InvalidPaths);
        }
    }
    for index in 0..second.frame_count() {
        let frame = second
            .frame(index)
            .ok_or(PathSimilarityError::InvalidPaths)?;
        if frame.len() != atom_count || frame.iter().flatten().any(|value| !value.is_finite()) {
            return Err(PathSimilarityError::InvalidPaths);
        }
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

#[cfg(test)]
#[path = "path_similarity_tests.rs"]
mod tests;
