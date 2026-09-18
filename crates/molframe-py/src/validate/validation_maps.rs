//! Native density-map validation and governed correlation projections.

use super::PyRealSpaceCorrelation;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::xtal_maps::{PyDensityMap, PyMapBoundary};
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn value_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}

fn map_boundary(value: PyMapBoundary) -> molframe::xtal::MapBoundary {
    match value {
        PyMapBoundary::Missing => molframe::xtal::MapBoundary::Missing,
        PyMapBoundary::Periodic => molframe::xtal::MapBoundary::Periodic,
    }
}

fn coordinate_values<'a>(values: &'a PyReadonlyArray2<'_, f64>) -> PyResult<&'a [[f64; 3]]> {
    let shape = values.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("positions must have shape (n, 3)"));
    }
    let values = values.as_slice().map_err(|_| {
        PyValueError::new_err(
            "positions must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    let (positions, remainder) = values.as_chunks::<3>();
    if remainder.is_empty() {
        Ok(positions)
    } else {
        Err(PyValueError::new_err("positions must have shape (n, 3)"))
    }
}

fn mask_values<'py>(mask: &'py PyReadonlyArray1<'py, bool>) -> PyResult<&'py [bool]> {
    mask.as_slice().map_err(|_| {
        PyValueError::new_err("mask must be C-contiguous; call numpy.ascontiguousarray explicitly")
    })
}

#[pyfunction]
pub(crate) fn real_space_map_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
) -> PyResult<PyRealSpaceCorrelation> {
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    py.detach(move || molframe::validate::real_space_map_correlation(&observed, &calculated))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn masked_real_space_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
    mask: PyReadonlyArray1<'_, bool>,
) -> PyResult<PyRealSpaceCorrelation> {
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    let mask = mask_values(&mask)?;
    py.detach(move || {
        molframe::validate::masked_real_space_correlation(&observed, &calculated, mask)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn sampled_real_space_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
    positions: PyReadonlyArray2<'_, f64>,
    boundary: PyMapBoundary,
) -> PyResult<PyRealSpaceCorrelation> {
    let positions = coordinate_values(&positions)?;
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    py.detach(move || {
        molframe::validate::sampled_real_space_correlation(
            &observed,
            &calculated,
            positions,
            map_boundary(boundary),
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

fn correlation_value(
    py: Python<'_>,
    value: molframe::validate::RealSpaceCorrelation,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyRealSpaceCorrelation::from(value))?.into_any())
}

#[pyfunction]
pub(crate) fn governed_real_space_map_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
    policy: &crate::query::PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::validate::governed_real_space_map_correlation(&observed, &calculated, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, correlation_value)
}

#[pyfunction]
pub(crate) fn governed_masked_real_space_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
    mask: PyReadonlyArray1<'_, bool>,
    policy: &crate::query::PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    let policy = policy.inner.clone();
    let mask = mask_values(&mask)?;
    let analysis = py
        .detach(move || {
            molframe::validate::governed_masked_real_space_correlation(
                &observed,
                &calculated,
                mask,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, correlation_value)
}

#[pyfunction]
pub(crate) fn governed_sampled_real_space_correlation(
    py: Python<'_>,
    observed: &PyDensityMap,
    calculated: &PyDensityMap,
    positions: PyReadonlyArray2<'_, f64>,
    policy: &crate::query::PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let positions = coordinate_values(&positions)?;
    let observed = observed.0.clone();
    let calculated = calculated.0.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::validate::governed_sampled_real_space_correlation(
                &observed,
                &calculated,
                positions,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, correlation_value)
}
