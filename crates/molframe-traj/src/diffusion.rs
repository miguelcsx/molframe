//! Diffusion-map embeddings over an explicit ensemble distance matrix.

use nalgebra::{DMatrix, SymmetricEigen};

use crate::ensemble::{EnsembleDistanceMatrix, EnsembleGeometryError};

/// A diffusion embedding with graph-connectivity diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionMap {
    /// Non-constant eigenvalues in descending order.
    pub eigenvalues: Box<[f64]>,
    /// Observation coordinates after applying diffusion time.
    pub coordinates: Box<[Box<[f64]>]>,
    /// Connected component label for every observation.
    pub graph_components: Box<[usize]>,
}

/// Builds an `exp(-d²/epsilon)` diffusion map with row normalization.
///
/// The constant first eigenvector is excluded. Returned coordinates describe
/// the supplied affinity geometry and are not asserted to be reaction coordinates.
///
/// # Errors
///
/// Returns an error for an invalid distance matrix, epsilon, or output dimension.
pub fn diffusion_map(
    distances: &EnsembleDistanceMatrix,
    epsilon: f64,
    diffusion_time: u32,
    dimensions: usize,
) -> Result<DiffusionMap, EnsembleGeometryError> {
    if distances.size < 2
        || distances.values.len() != distances.size * distances.size
        || dimensions == 0
    {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    if !epsilon.is_finite() || epsilon <= 0.0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let size = distances.size;
    let mut affinity = DMatrix::zeros(size, size);
    for row in 0..size {
        for column in 0..size {
            let distance = distances.values[row * size + column];
            if !distance.is_finite() || distance < 0.0 {
                return Err(EnsembleGeometryError::InvalidParameter);
            }
            affinity[(row, column)] = (-(distance * distance) / epsilon).exp();
        }
    }
    let degrees: Vec<f64> = (0..size)
        .map(|row| (0..size).map(|column| affinity[(row, column)]).sum())
        .collect();
    let symmetric = DMatrix::from_fn(size, size, |row, column| {
        affinity[(row, column)] / (degrees[row] * degrees[column]).sqrt()
    });
    let decomposition = SymmetricEigen::new(symmetric);
    let mut indices: Vec<usize> = (0..size).collect();
    indices.sort_by(|&left, &right| {
        decomposition.eigenvalues[right]
            .total_cmp(&decomposition.eigenvalues[left])
            .then(left.cmp(&right))
    });
    let selected: Vec<usize> = indices.into_iter().skip(1).take(dimensions).collect();
    let eigenvalues: Vec<f64> = selected
        .iter()
        .map(|&index| decomposition.eigenvalues[index])
        .collect();
    let coordinates = (0..size)
        .map(|row| {
            selected
                .iter()
                .map(|&index| {
                    let right_eigenvector =
                        decomposition.eigenvectors[(row, index)] / degrees[row].sqrt();
                    right_eigenvector
                        * decomposition.eigenvalues[index].powf(f64::from(diffusion_time))
                })
                .collect::<Vec<f64>>()
                .into_boxed_slice()
        })
        .collect();
    Ok(DiffusionMap {
        eigenvalues: eigenvalues.into_boxed_slice(),
        coordinates,
        graph_components: connected_components(&affinity).into_boxed_slice(),
    })
}

fn connected_components(affinity: &DMatrix<f64>) -> Vec<usize> {
    let size = affinity.nrows();
    let mut labels = vec![usize::MAX; size];
    let mut label = 0;
    for start in 0..size {
        if labels[start] != usize::MAX {
            continue;
        }
        labels[start] = label;
        let mut stack = vec![start];
        while let Some(row) = stack.pop() {
            for column in 0..size {
                if row != column && affinity[(row, column)] > 0.0 && labels[column] == usize::MAX {
                    labels[column] = label;
                    stack.push(column);
                }
            }
        }
        label += 1;
    }
    labels
}

#[cfg(test)]
#[path = "diffusion_tests.rs"]
mod tests;
