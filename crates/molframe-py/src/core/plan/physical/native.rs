//! Lower atom-aligned physical inputs into borrowed facade slots.

use super::PyPhysicalOperation;
use numpy::PyReadonlyArray1;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;
use std::collections::HashMap;

pub(crate) fn to_native<'py>(
    py: Python<'py>,
    operation: &PyPhysicalOperation,
    scalar_arrays: &mut Vec<PyReadonlyArray1<'py, f64>>,
    float_arrays: &mut Vec<PyReadonlyArray1<'py, f32>>,
    scalar_slots: &mut HashMap<usize, usize>,
    float_slots: &mut HashMap<usize, usize>,
) -> PyResult<molframe::PhysicalRequest> {
    Ok(match operation {
        PyPhysicalOperation::RadialDistribution {
            left,
            right,
            options,
            periodic,
            policy,
        } => molframe::PhysicalRequest::RadialDistribution {
            left: left.clone(),
            right: right.clone(),
            options: *options,
            periodic: *periodic,
            policy: policy.clone(),
        },
        PyPhysicalOperation::CoordinationNumbers {
            left,
            right,
            minimum_distance,
            maximum_distance,
            backend,
            periodic,
            policy,
        } => molframe::PhysicalRequest::CoordinationNumbers {
            left: left.clone(),
            right: right.clone(),
            minimum_distance: *minimum_distance,
            maximum_distance: *maximum_distance,
            backend: *backend,
            periodic: *periodic,
            policy: policy.clone(),
        },
        PyPhysicalOperation::Leaflets {
            sites,
            options,
            periodic,
            policy,
        } => molframe::PhysicalRequest::Leaflets {
            sites: sites.clone(),
            options: *options,
            periodic: *periodic,
            policy: policy.clone(),
        },
        PyPhysicalOperation::LinearDensity {
            weights,
            options,
            policy,
        } => molframe::PhysicalRequest::LinearDensity {
            weights: retain_scalar(py, weights, scalar_arrays, scalar_slots)?,
            options: *options,
            policy: policy.clone(),
        },
        PyPhysicalOperation::DensityMap {
            weights,
            spec,
            policy,
        } => molframe::PhysicalRequest::DensityMap {
            weights: retain_scalar(py, weights, scalar_arrays, scalar_slots)?,
            spec: *spec,
            policy: policy.clone(),
        },
        PyPhysicalOperation::PoreProfile {
            radii,
            options,
            policy,
        } => molframe::PhysicalRequest::PoreProfile {
            radii: retain_float(py, radii, float_arrays, float_slots)?,
            options: *options,
            policy: policy.clone(),
        },
        PyPhysicalOperation::SurfaceContacts {
            radii,
            tolerance,
            probe,
            density,
            minimum_area,
            backend,
            policy,
        } => molframe::PhysicalRequest::SurfaceContacts {
            radii: retain_float(py, radii, float_arrays, float_slots)?,
            tolerance: *tolerance,
            probe: *probe,
            density: *density,
            minimum_area: *minimum_area,
            backend: *backend,
            policy: policy.clone(),
        },
    })
}

fn retain_scalar<'py>(
    py: Python<'py>,
    value: &Py<PyAny>,
    arrays: &mut Vec<PyReadonlyArray1<'py, f64>>,
    slots: &mut HashMap<usize, usize>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray1<'py, f64>>()?;
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "physical scalar inputs must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let slot = arrays.len();
    arrays.push(array);
    slots.insert(key, slot);
    Ok(slot)
}

fn retain_float<'py>(
    py: Python<'py>,
    value: &Py<PyAny>,
    arrays: &mut Vec<PyReadonlyArray1<'py, f32>>,
    slots: &mut HashMap<usize, usize>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray1<'py, f32>>()?;
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "physical float inputs must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let slot = arrays.len();
    arrays.push(array);
    slots.insert(key, slot);
    Ok(slot)
}
