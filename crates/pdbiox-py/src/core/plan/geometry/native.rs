//! Lower typed geometry nodes into borrowed native inputs.

use super::super::dispatch::require_operation_tag;
use super::model::PyGeometryOperation;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::HashMap;

pub(crate) fn serialized_request(
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PyGeometryOperation>> {
    let Some(value) = config.get_item("operation")? else {
        return Ok(None);
    };
    let operation: &str = value.extract()?;
    let request = match operation {
        "centroid" => position_from_dict(config, operation, |positions| {
            PyGeometryOperation::Centroid { positions }
        })?,
        "distance_matrix" => position_from_dict(config, operation, |positions| {
            PyGeometryOperation::DistanceMatrix { positions }
        })?,
        "centre_of_mass" => weighted_from_dict(config, operation, |positions, masses| {
            PyGeometryOperation::CentreOfMass { positions, masses }
        })?,
        "radius_of_gyration" => weighted_from_dict(config, operation, |positions, masses| {
            PyGeometryOperation::RadiusOfGyration { positions, masses }
        })?,
        "inertia_tensor" => weighted_from_dict(config, operation, |positions, masses| {
            PyGeometryOperation::InertiaTensor { positions, masses }
        })?,
        "principal_axes" => principal_from_dict(config)?,
        "asphericity" => eigen_from_dict(config, operation, |positions, options| {
            PyGeometryOperation::Asphericity { positions, options }
        })?,
        "gyration_axes" => eigen_from_dict(config, operation, |positions, options| {
            PyGeometryOperation::GyrationAxes { positions, options }
        })?,
        "distance_matrix_between" => {
            let request = PyGeometryOperation::DistanceMatrixBetween {
                left: super::super::required(config, "left")?.unbind(),
                right: super::super::required(config, "right")?.unbind(),
            };
            validate_operation(config.py(), &request)?;
            request
        }
        "rmsf" => {
            let request = PyGeometryOperation::Rmsf {
                frames: super::super::required(config, "frames")?.unbind(),
            };
            validate_operation(config.py(), &request)?;
            request
        }
        _ => return Ok(None),
    };
    Ok(Some(request))
}

pub(crate) fn typed_request(value: &Bound<'_, PyAny>) -> PyResult<Option<PyGeometryOperation>> {
    super::classes::typed_classes(value)
}

pub(crate) fn to_native<'py>(
    py: Python<'py>,
    operation: &PyGeometryOperation,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    scalars: &mut Vec<PyReadonlyArray1<'py, f64>>,
    frames: &mut Vec<PyReadonlyArray3<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    frame_slots: &mut HashMap<usize, usize>,
) -> PyResult<pdbiox::GeometryRequest> {
    let mut positions = |value: &Py<PyAny>| retain_positions(py, value, arrays, coordinate_slots);
    let mut masses = |value: &Option<Py<PyAny>>| retain_masses(py, value.as_ref(), scalars);
    Ok(match operation {
        PyGeometryOperation::Centroid { positions: value } => pdbiox::GeometryRequest::Centroid {
            positions: positions(value)?,
        },
        PyGeometryOperation::CentreOfMass {
            positions: value,
            masses: weights,
        } => pdbiox::GeometryRequest::CentreOfMass {
            positions: positions(value)?,
            masses: masses(weights)?,
        },
        PyGeometryOperation::RadiusOfGyration {
            positions: value,
            masses: weights,
        } => pdbiox::GeometryRequest::RadiusOfGyration {
            positions: positions(value)?,
            masses: masses(weights)?,
        },
        PyGeometryOperation::InertiaTensor {
            positions: value,
            masses: weights,
        } => pdbiox::GeometryRequest::InertiaTensor {
            positions: positions(value)?,
            masses: masses(weights)?,
        },
        PyGeometryOperation::PrincipalAxes {
            positions: value,
            masses: weights,
            options,
        } => pdbiox::GeometryRequest::PrincipalAxes {
            positions: positions(value)?,
            masses: masses(weights)?,
            options: *options,
        },
        PyGeometryOperation::Asphericity {
            positions: value,
            options,
        } => pdbiox::GeometryRequest::Asphericity {
            positions: positions(value)?,
            options: *options,
        },
        PyGeometryOperation::GyrationAxes {
            positions: value,
            options,
        } => pdbiox::GeometryRequest::GyrationAxes {
            positions: positions(value)?,
            options: *options,
        },
        PyGeometryOperation::DistanceMatrix { positions: value } => {
            pdbiox::GeometryRequest::DistanceMatrix {
                positions: positions(value)?,
            }
        }
        PyGeometryOperation::DistanceMatrixBetween { left, right } => {
            pdbiox::GeometryRequest::DistanceMatrixBetween {
                left: positions(left)?,
                right: positions(right)?,
            }
        }
        PyGeometryOperation::Rmsf { frames: value } => pdbiox::GeometryRequest::Rmsf {
            frames: super::super::lowering::retain_frame(py, frames, frame_slots, value)?,
        },
    })
}

/// Validates the fixed-shape, zero-copy input contract of a geometry node.
///
/// This only inspects dtype, shape, contiguity, and aligned lengths. Coordinate
/// values remain borrowed and are read by the native kernel at execution time.
pub(crate) fn validate_operation(py: Python<'_>, operation: &PyGeometryOperation) -> PyResult<()> {
    match operation {
        PyGeometryOperation::Centroid { positions }
        | PyGeometryOperation::DistanceMatrix { positions }
        | PyGeometryOperation::Asphericity { positions, .. }
        | PyGeometryOperation::GyrationAxes { positions, .. } => {
            validate_positions(py, positions)?;
        }
        PyGeometryOperation::CentreOfMass { positions, masses }
        | PyGeometryOperation::RadiusOfGyration { positions, masses }
        | PyGeometryOperation::InertiaTensor { positions, masses }
        | PyGeometryOperation::PrincipalAxes {
            positions, masses, ..
        } => {
            let count = validate_positions(py, positions)?;
            validate_masses(py, masses.as_ref(), count)?;
        }
        PyGeometryOperation::DistanceMatrixBetween { left, right } => {
            validate_positions(py, left)?;
            validate_positions(py, right)?;
        }
        PyGeometryOperation::Rmsf { frames } => validate_frames(py, frames)?,
    }
    Ok(())
}

fn validate_positions(py: Python<'_>, value: &Py<PyAny>) -> PyResult<usize> {
    let array = value.bind(py).extract::<PyReadonlyArray2<'_, f32>>()?;
    Ok(crate::geometry::borrowed_coordinates(&array)?.len())
}

fn validate_masses(
    py: Python<'_>,
    value: Option<&Py<PyAny>>,
    coordinate_count: usize,
) -> PyResult<()> {
    let Some(value) = value else {
        return Ok(());
    };
    let array = value.bind(py).extract::<PyReadonlyArray1<'_, f64>>()?;
    let masses = array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "masses must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    if masses.len() == coordinate_count {
        Ok(())
    } else {
        Err(PyValueError::new_err(
            "masses must contain one value per coordinate",
        ))
    }
}

fn validate_frames(py: Python<'_>, value: &Py<PyAny>) -> PyResult<()> {
    let frames = value.bind(py).extract::<PyReadonlyArray3<'_, f32>>()?;
    crate::intrinsic::borrowed_frame_input(&frames).map(|_| ())
}

fn retain_positions<'py>(
    py: Python<'py>,
    value: &Py<PyAny>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = coordinate_slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray2<'py, f32>>()?;
    crate::geometry::borrowed_coordinates(&array)?;
    let slot = arrays.len();
    arrays.push(array);
    coordinate_slots.insert(key, slot);
    Ok(slot)
}

fn retain_masses<'py>(
    py: Python<'py>,
    value: Option<&Py<PyAny>>,
    scalars: &mut Vec<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Option<usize>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let array = value.bind(py).extract::<PyReadonlyArray1<'py, f64>>()?;
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "masses must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let slot = scalars.len();
    scalars.push(array);
    Ok(Some(slot))
}

pub(crate) fn position_from_dict(
    config: &Bound<'_, PyDict>,
    tag: &str,
    make: fn(Py<PyAny>) -> PyGeometryOperation,
) -> PyResult<PyGeometryOperation> {
    require_operation_tag(config, tag)?;
    let operation = make(super::super::required(config, "positions")?.unbind());
    validate_operation(config.py(), &operation)?;
    Ok(operation)
}

pub(crate) fn weighted_from_dict(
    config: &Bound<'_, PyDict>,
    tag: &str,
    make: fn(Py<PyAny>, Option<Py<PyAny>>) -> PyGeometryOperation,
) -> PyResult<PyGeometryOperation> {
    require_operation_tag(config, tag)?;
    let masses = match config.get_item("masses")? {
        Some(value) if !value.is_none() => Some(value.unbind()),
        _ => None,
    };
    let operation = make(
        super::super::required(config, "positions")?.unbind(),
        masses,
    );
    validate_operation(config.py(), &operation)?;
    Ok(operation)
}

pub(crate) fn eigen_from_dict(
    config: &Bound<'_, PyDict>,
    tag: &str,
    make: fn(Py<PyAny>, pdbiox::EigenOptions) -> PyGeometryOperation,
) -> PyResult<PyGeometryOperation> {
    require_operation_tag(config, tag)?;
    let operation = make(
        super::super::required(config, "positions")?.unbind(),
        eigen_options(config)?,
    );
    validate_operation(config.py(), &operation)?;
    Ok(operation)
}

pub(crate) fn principal_from_dict(config: &Bound<'_, PyDict>) -> PyResult<PyGeometryOperation> {
    require_operation_tag(config, "principal_axes")?;
    let masses = match config.get_item("masses")? {
        Some(value) if !value.is_none() => Some(value.unbind()),
        _ => None,
    };
    let operation = PyGeometryOperation::PrincipalAxes {
        positions: super::super::required(config, "positions")?.unbind(),
        masses,
        options: eigen_options(config)?,
    };
    validate_operation(config.py(), &operation)?;
    Ok(operation)
}

fn eigen_options(config: &Bound<'_, PyDict>) -> PyResult<pdbiox::EigenOptions> {
    match config.get_item("options")? {
        Some(value) => Ok(value.extract::<crate::geometry::PyEigenOptions>()?.inner),
        None => Ok(pdbiox::EigenOptions::standard()),
    }
}
