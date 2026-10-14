//! Allocation ceilings for PCA requests, checked before any buffer exists.

use crate::ensemble::EnsembleGeometryError;
use crate::pca::CartesianFit;
use std::mem::size_of;

pub(crate) fn check_pca_request(
    observations: usize,
    features: usize,
    component_count: usize,
    memory_limit: usize,
    fit_memory: FitMemory,
) -> Result<(), EnsembleGeometryError> {
    if observations < 2 || features == 0 || component_count == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let components = component_count.min(observations - 1).min(features);
    let square = observations.min(features);
    let data = observations
        .checked_mul(features)
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let eigensolver = square
        .checked_mul(square)
        .and_then(|entries| entries.checked_mul(2))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let result = features
        .checked_mul(components)
        .and_then(|component_entries| {
            observations
                .checked_mul(components)
                .and_then(|projection_entries| component_entries.checked_add(projection_entries))
        })
        .and_then(|entries| entries.checked_add(features))
        .and_then(|entries| entries.checked_add(components))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let pca_bytes = data
        .checked_add(eigensolver)
        .and_then(|entries| entries.checked_add(result))
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let matrix_bytes = data
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(fit_memory.retained_bytes))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let required = pca_bytes
        .max(matrix_bytes)
        .max(fit_memory.construction_peak);
    if required > memory_limit {
        return Err(EnsembleGeometryError::MemoryLimit);
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(crate) struct FitMemory {
    pub(crate) retained_bytes: usize,
    pub(crate) construction_peak: usize,
}

impl FitMemory {
    pub(crate) const NONE: Self = Self {
        retained_bytes: 0,
        construction_peak: 0,
    };
}

pub(crate) fn fit_memory(
    atoms: usize,
    fit: CartesianFit<'_>,
) -> Result<FitMemory, EnsembleGeometryError> {
    let CartesianFit::IterativeMean {
        tolerance,
        max_iterations,
    } = fit
    else {
        return Ok(FitMemory::NONE);
    };
    if !tolerance.is_finite() || tolerance <= 0.0 || max_iterations == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let retained_bytes = atoms
        .checked_mul(size_of::<[f32; 3]>())
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let accumulator_bytes = atoms
        .checked_mul(size_of::<[f64; 3]>())
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let construction_peak = retained_bytes
        .checked_mul(2)
        .and_then(|bytes| bytes.checked_add(accumulator_bytes))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    Ok(FitMemory {
        retained_bytes,
        construction_peak,
    })
}
