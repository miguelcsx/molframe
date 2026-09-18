//! Declarative surface operations over borrowed `NumPy` buffers.

use super::dispatch::require_operation_tag;
use super::lowering::retain_coordinate;
use super::{Operation, PyPlan, execute_native};
use crate::geometry::borrowed_coordinates;
use crate::surface_types::PyBuriedSurface;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::HashMap;

/// Owned Python configuration for one native surface request.
pub(crate) enum PySurfaceOperation {
    SolventAccessibleSurface {
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    },
    BuriedSurface {
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        first: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    },
}

impl PySurfaceOperation {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::SolventAccessibleSurface {
                positions,
                radii,
                probe,
                sample_points,
            } => Self::SolventAccessibleSurface {
                positions: positions.clone_ref(py),
                radii: radii.clone_ref(py),
                probe: *probe,
                sample_points: *sample_points,
            },
            Self::BuriedSurface {
                positions,
                radii,
                first,
                probe,
                sample_points,
            } => Self::BuriedSurface {
                positions: positions.clone_ref(py),
                radii: radii.clone_ref(py),
                first: first.clone_ref(py),
                probe: *probe,
                sample_points: *sample_points,
            },
        }
    }
}

#[pyclass(name = "Sasa", frozen)]
pub(crate) struct PySasa {
    pub(crate) operation: PySurfaceOperation,
}

#[pymethods]
impl PySasa {
    #[new]
    #[pyo3(signature = (*, positions, radii, probe, sample_points))]
    fn new(
        py: Python<'_>,
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    ) -> PyResult<Self> {
        Self::build(py, positions, radii, probe, sample_points)
    }

    fn execute<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray1<f64>>> {
        match execute_operation(py, &self.operation)? {
            molframe::SurfaceValue::SolventAccessibleSurface(values) => {
                use numpy::IntoPyArray;
                Ok(numpy::ndarray::Array1::from_vec(values).into_pyarray(py))
            }
            molframe::SurfaceValue::BuriedSurface(_) => Err(incompatible_result()),
        }
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "sasa")?;
        Self::build(
            py,
            required(config, "positions")?.unbind(),
            required(config, "radii")?.unbind(),
            required(config, "probe")?.extract()?,
            required(config, "sample_points")?.extract()?,
        )
    }

    fn __repr__(&self) -> String {
        "Sasa(positions=<array>, radii=<array>, probe=..., sample_points=...)".to_owned()
    }
}

impl PySasa {
    fn build(
        py: Python<'_>,
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    ) -> PyResult<Self> {
        validate_inputs(py, &positions, &radii, None, probe, sample_points)?;
        Ok(Self {
            operation: PySurfaceOperation::SolventAccessibleSurface {
                positions,
                radii,
                probe,
                sample_points,
            },
        })
    }
}

#[pyclass(name = "BuriedSurfaceOp", frozen)]
pub(crate) struct PyBuriedSurfaceOperation {
    pub(crate) operation: PySurfaceOperation,
}

#[pymethods]
impl PyBuriedSurfaceOperation {
    #[new]
    #[pyo3(signature = (*, positions, radii, first, probe, sample_points))]
    fn new(
        py: Python<'_>,
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        first: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    ) -> PyResult<Self> {
        Self::build(py, positions, radii, first, probe, sample_points)
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyBuriedSurface> {
        match execute_operation(py, &self.operation)? {
            molframe::SurfaceValue::BuriedSurface(value) => Ok(value.into()),
            molframe::SurfaceValue::SolventAccessibleSurface(_) => Err(incompatible_result()),
        }
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "buried_surface")?;
        Self::build(
            py,
            required(config, "positions")?.unbind(),
            required(config, "radii")?.unbind(),
            required(config, "first")?.unbind(),
            required(config, "probe")?.extract()?,
            required(config, "sample_points")?.extract()?,
        )
    }

    fn __repr__(&self) -> String {
        "BuriedSurfaceOp(positions=<array>, radii=<array>, first=<mask>, probe=..., sample_points=...)"
            .to_owned()
    }
}

impl PyBuriedSurfaceOperation {
    fn build(
        py: Python<'_>,
        positions: Py<PyAny>,
        radii: Py<PyAny>,
        first: Py<PyAny>,
        probe: f32,
        sample_points: u16,
    ) -> PyResult<Self> {
        validate_inputs(py, &positions, &radii, Some(&first), probe, sample_points)?;
        Ok(Self {
            operation: PySurfaceOperation::BuriedSurface {
                positions,
                radii,
                first,
                probe,
                sample_points,
            },
        })
    }
}

fn execute_operation(
    py: Python<'_>,
    operation: &PySurfaceOperation,
) -> PyResult<molframe::SurfaceValue> {
    let plan = PyPlan {
        operations: vec![(
            "result".to_owned(),
            Operation::Surface(operation.clone_ref(py)),
        )],
    };
    let result = execute_native(&plan, py, None, None)?;
    let Some(entry) = result.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native surface plan returned no result",
        ));
    };
    match entry.value {
        molframe::PlanValue::Surface(value) => Ok(*value),
        _ => Err(incompatible_result()),
    }
}

pub(crate) fn to_native<'py>(
    py: Python<'py>,
    operation: &PySurfaceOperation,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    floats: &mut Vec<PyReadonlyArray1<'py, f32>>,
    masks: &mut Vec<PyReadonlyArray1<'py, bool>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    float_slots: &mut HashMap<usize, usize>,
    mask_slots: &mut HashMap<usize, usize>,
) -> PyResult<molframe::SurfaceRequest> {
    match operation {
        PySurfaceOperation::SolventAccessibleSurface {
            positions,
            radii,
            probe,
            sample_points,
        } => Ok(molframe::SurfaceRequest::SolventAccessibleSurface {
            positions: retain_coordinate(py, arrays, coordinate_slots, positions)?,
            radii: retain_float(py, radii, floats, float_slots)?,
            probe: *probe,
            sample_points: *sample_points,
        }),
        PySurfaceOperation::BuriedSurface {
            positions,
            radii,
            first,
            probe,
            sample_points,
        } => Ok(molframe::SurfaceRequest::BuriedSurface {
            positions: retain_coordinate(py, arrays, coordinate_slots, positions)?,
            radii: retain_float(py, radii, floats, float_slots)?,
            first: retain_mask(py, first, masks, mask_slots)?,
            probe: *probe,
            sample_points: *sample_points,
        }),
    }
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
    contiguous_float(&array)?;
    let slot = arrays.len();
    arrays.push(array);
    slots.insert(key, slot);
    Ok(slot)
}

fn retain_mask<'py>(
    py: Python<'py>,
    value: &Py<PyAny>,
    arrays: &mut Vec<PyReadonlyArray1<'py, bool>>,
    slots: &mut HashMap<usize, usize>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray1<'py, bool>>()?;
    contiguous_mask(&array)?;
    let slot = arrays.len();
    arrays.push(array);
    slots.insert(key, slot);
    Ok(slot)
}

fn validate_inputs(
    py: Python<'_>,
    positions: &Py<PyAny>,
    radii: &Py<PyAny>,
    first: Option<&Py<PyAny>>,
    probe: f32,
    sample_points: u16,
) -> PyResult<()> {
    let coordinates = positions.bind(py).extract::<PyReadonlyArray2<'_, f32>>()?;
    let coordinate_count = borrowed_coordinates(&coordinates)?.len();
    let radii = radii.bind(py).extract::<PyReadonlyArray1<'_, f32>>()?;
    contiguous_float(&radii)?;
    if radii.len()? != coordinate_count {
        return Err(PyValueError::new_err(
            "radii must contain one value per coordinate",
        ));
    }
    if let Some(first) = first {
        let first = first.bind(py).extract::<PyReadonlyArray1<'_, bool>>()?;
        contiguous_mask(&first)?;
        if first.len()? != coordinate_count {
            return Err(PyValueError::new_err(
                "first must contain one mask value per coordinate",
            ));
        }
    }
    if !probe.is_finite() || probe < 0.0 {
        return Err(PyValueError::new_err(
            "probe must be finite and non-negative",
        ));
    }
    if sample_points == 0 {
        return Err(PyValueError::new_err("sample_points must be positive"));
    }
    Ok(())
}

fn contiguous_float(array: &PyReadonlyArray1<'_, f32>) -> PyResult<()> {
    array.as_slice().map(|_| ()).map_err(|_| {
        PyValueError::new_err(
            "surface radii must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}

fn contiguous_mask(array: &PyReadonlyArray1<'_, bool>) -> PyResult<()> {
    array.as_slice().map(|_| ()).map_err(|_| {
        PyValueError::new_err(
            "surface masks must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}

pub(crate) fn typed_request(value: &Bound<'_, PyAny>) -> Option<PySurfaceOperation> {
    if let Ok(value) = value.extract::<PyRef<'_, PySasa>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyBuriedSurfaceOperation>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    None
}

pub(crate) fn serialized_request(
    py: Python<'_>,
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PySurfaceOperation>> {
    let Some(operation) = config.get_item("operation")? else {
        return Ok(None);
    };
    match operation.extract::<&str>()? {
        "sasa" => Ok(Some(PySasa::from_dict(py, config)?.operation)),
        "buried_surface" => Ok(Some(
            PyBuriedSurfaceOperation::from_dict(py, config)?.operation,
        )),
        _ => Ok(None),
    }
}

pub(super) fn explain_operation<'py>(
    py: Python<'py>,
    operation: &PySurfaceOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation_tag(operation))?;
    result.set_item("execution", "native")?;
    result.set_item("gil_released", true)?;
    result.set_item("requires_c_contiguous", true)?;
    result.set_item("materializes", materialization(operation))?;
    result.set_item("complexity", "O(atoms × samples × local density)")?;
    Ok(result)
}

pub(super) fn operation_to_dict<'py>(
    py: Python<'py>,
    operation: &PySurfaceOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    match operation {
        PySurfaceOperation::SolventAccessibleSurface {
            positions,
            radii,
            probe,
            sample_points,
        } => {
            result.set_item("operation", "sasa")?;
            result.set_item("positions", positions.bind(py))?;
            result.set_item("radii", radii.bind(py))?;
            result.set_item("probe", probe)?;
            result.set_item("sample_points", sample_points)?;
        }
        PySurfaceOperation::BuriedSurface {
            positions,
            radii,
            first,
            probe,
            sample_points,
        } => {
            result.set_item("operation", "buried_surface")?;
            result.set_item("positions", positions.bind(py))?;
            result.set_item("radii", radii.bind(py))?;
            result.set_item("first", first.bind(py))?;
            result.set_item("probe", probe)?;
            result.set_item("sample_points", sample_points)?;
        }
    }
    Ok(result)
}

fn operation_tag(operation: &PySurfaceOperation) -> &'static str {
    match operation {
        PySurfaceOperation::SolventAccessibleSurface { .. } => "sasa",
        PySurfaceOperation::BuriedSurface { .. } => "buried_surface",
    }
}

fn materialization(operation: &PySurfaceOperation) -> &'static str {
    match operation {
        PySurfaceOperation::SolventAccessibleSurface { .. } => "per-atom area array",
        PySurfaceOperation::BuriedSurface { .. } => "surface area summary",
    }
}

fn required<'py>(config: &'py Bound<'py, PyDict>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    config
        .get_item(name)?
        .ok_or_else(|| PyValueError::new_err(format!("missing operation field: {name}")))
}

fn incompatible_result() -> PyErr {
    PyTypeError::new_err("native surface plan returned an incompatible result")
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySasa>()?;
    module.add_class::<PyBuriedSurfaceOperation>()?;
    Ok(())
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
