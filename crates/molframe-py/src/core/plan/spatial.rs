//! Declarative fixed-radius spatial operations over borrowed `NumPy` arrays.

use super::dispatch::require_operation_tag;
use super::{Operation, PyPlan, execute_native};
use crate::spatial::{
    PyNeighborTable, PySpatialSearchOptions, spatial_arrays::SpatialBindingError,
};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::HashMap;

#[path = "spatial_helpers.rs"]
mod helpers;
use helpers::{
    incompatible_result, optional_object, owned_selection, periodic_box, position_count, required,
    retain_positions, set_optional_object, spatial_options, validate_cell, validate_cutoff,
    validate_indices, validate_optional_indices,
};

/// Owned Python configuration for one native spatial request.
pub(crate) enum PySpatialOperation {
    NeighborPairs {
        positions: Py<PyAny>,
        left: Option<Py<PyAny>>,
        right: Option<Py<PyAny>>,
        cutoff: f32,
        options: PySpatialSearchOptions,
        cell: Option<Py<PyAny>>,
    },
    AtomsWithin {
        positions: Py<PyAny>,
        targets: Py<PyAny>,
        query: Option<Py<PyAny>>,
        cutoff: f32,
        options: PySpatialSearchOptions,
        cell: Option<Py<PyAny>>,
    },
}

impl PySpatialOperation {
    fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::NeighborPairs {
                positions,
                left,
                right,
                cutoff,
                options,
                cell,
            } => Self::NeighborPairs {
                positions: positions.clone_ref(py),
                left: left.as_ref().map(|value| value.clone_ref(py)),
                right: right.as_ref().map(|value| value.clone_ref(py)),
                cutoff: *cutoff,
                options: *options,
                cell: cell.as_ref().map(|value| value.clone_ref(py)),
            },
            Self::AtomsWithin {
                positions,
                targets,
                query,
                cutoff,
                options,
                cell,
            } => Self::AtomsWithin {
                positions: positions.clone_ref(py),
                targets: targets.clone_ref(py),
                query: query.as_ref().map(|value| value.clone_ref(py)),
                cutoff: *cutoff,
                options: *options,
                cell: cell.as_ref().map(|value| value.clone_ref(py)),
            },
        }
    }
}

#[pyclass(name = "NeighborPairs", frozen)]
pub(crate) struct PyNeighborPairs {
    pub(crate) operation: PySpatialOperation,
}

#[pymethods]
impl PyNeighborPairs {
    #[new]
    #[pyo3(signature = (*, positions, cutoff, left=None, right=None, options=None, cell=None))]
    fn new(
        py: Python<'_>,
        positions: Py<PyAny>,
        cutoff: f32,
        left: Option<Py<PyAny>>,
        right: Option<Py<PyAny>>,
        options: Option<Py<PyAny>>,
        cell: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        Self::build(py, positions, left, right, cutoff, options, cell)
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyNeighborTable> {
        match execute_operation(py, &self.operation)? {
            molframe::SpatialValue::NeighborPairs(values) => neighbor_table_from_pairs(py, values),
            molframe::SpatialValue::AtomsWithin(_) => Err(incompatible_result()),
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
        require_operation_tag(config, "neighbor_pairs")?;
        let positions = required(config, "positions")?.unbind();
        Self::build(
            py,
            positions,
            optional_object(config, "left")?,
            optional_object(config, "right")?,
            required(config, "cutoff")?.extract()?,
            optional_object(config, "options")?,
            optional_object(config, "cell")?,
        )
    }

    fn __repr__(&self) -> String {
        "NeighborPairs(positions=<array>, cutoff=..., options=...)".to_owned()
    }
}

impl PyNeighborPairs {
    fn build(
        py: Python<'_>,
        positions: Py<PyAny>,
        left: Option<Py<PyAny>>,
        right: Option<Py<PyAny>>,
        cutoff: f32,
        options: Option<Py<PyAny>>,
        cell: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        validate_cutoff(cutoff)?;
        let count = position_count(py, &positions)?;
        validate_optional_indices(py, left.as_ref(), count)?;
        validate_optional_indices(py, right.as_ref(), count)?;
        let options = spatial_options(py, options.as_ref())?;
        validate_cell(py, cell.as_ref())?;
        Ok(Self {
            operation: PySpatialOperation::NeighborPairs {
                positions,
                left,
                right,
                cutoff,
                options,
                cell,
            },
        })
    }
}

#[pyclass(name = "AtomsWithin", frozen)]
pub(crate) struct PyAtomsWithin {
    pub(crate) operation: PySpatialOperation,
}

#[pymethods]
impl PyAtomsWithin {
    #[new]
    #[pyo3(signature = (*, positions, targets, cutoff, query=None, options=None, cell=None))]
    fn new(
        py: Python<'_>,
        positions: Py<PyAny>,
        targets: Py<PyAny>,
        cutoff: f32,
        query: Option<Py<PyAny>>,
        options: Option<Py<PyAny>>,
        cell: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        Self::build(py, positions, targets, query, cutoff, options, cell)
    }

    fn execute<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u32>>> {
        match execute_operation(py, &self.operation)? {
            molframe::SpatialValue::AtomsWithin(value) => Ok(atom_indices_to_python(py, value)),
            molframe::SpatialValue::NeighborPairs(_) => Err(incompatible_result()),
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
        require_operation_tag(config, "atoms_within")?;
        let positions = required(config, "positions")?.unbind();
        let targets = required(config, "targets")?.unbind();
        Self::build(
            py,
            positions,
            targets,
            optional_object(config, "query")?,
            required(config, "cutoff")?.extract()?,
            optional_object(config, "options")?,
            optional_object(config, "cell")?,
        )
    }

    fn __repr__(&self) -> String {
        "AtomsWithin(positions=<array>, targets=..., cutoff=..., options=...)".to_owned()
    }
}

impl PyAtomsWithin {
    fn build(
        py: Python<'_>,
        positions: Py<PyAny>,
        targets: Py<PyAny>,
        query: Option<Py<PyAny>>,
        cutoff: f32,
        options: Option<Py<PyAny>>,
        cell: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        validate_cutoff(cutoff)?;
        let count = position_count(py, &positions)?;
        validate_indices(py, &targets, count)?;
        validate_optional_indices(py, query.as_ref(), count)?;
        let options = spatial_options(py, options.as_ref())?;
        validate_cell(py, cell.as_ref())?;
        Ok(Self {
            operation: PySpatialOperation::AtomsWithin {
                positions,
                targets,
                query,
                cutoff,
                options,
                cell,
            },
        })
    }
}

fn execute_operation(
    py: Python<'_>,
    operation: &PySpatialOperation,
) -> PyResult<molframe::SpatialValue> {
    let plan = PyPlan {
        operations: vec![(
            "result".to_owned(),
            Operation::Spatial(operation.clone_ref(py)),
        )],
    };
    let result = execute_native(&plan, py, None, None)?;
    let Some(entry) = result.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native spatial plan returned no result",
        ));
    };
    match entry.value {
        molframe::PlanValue::Spatial(value) => Ok(*value),
        _ => Err(incompatible_result()),
    }
}

pub(crate) fn to_native<'py>(
    py: Python<'py>,
    operation: &PySpatialOperation,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<molframe::SpatialRequest> {
    match operation {
        PySpatialOperation::NeighborPairs {
            positions,
            left,
            right,
            cutoff,
            options,
            cell,
        } => {
            let slot = retain_positions(py, positions, arrays, coordinate_slots)?;
            let count = arrays[slot].shape()[0];
            Ok(molframe::SpatialRequest::NeighborPairs {
                positions: slot,
                left: owned_selection(py, left.as_ref(), count)?,
                right: owned_selection(py, right.as_ref(), count)?,
                cutoff: *cutoff,
                options: options.inner(),
                periodic: periodic_box(py, cell.as_ref())?,
            })
        }
        PySpatialOperation::AtomsWithin {
            positions,
            targets,
            query,
            cutoff,
            options,
            cell,
        } => {
            let slot = retain_positions(py, positions, arrays, coordinate_slots)?;
            let count = arrays[slot].shape()[0];
            Ok(molframe::SpatialRequest::AtomsWithin {
                positions: slot,
                query: owned_selection(py, query.as_ref(), count)?,
                target: owned_selection(py, Some(targets), count)?,
                cutoff: *cutoff,
                options: options.inner(),
                periodic: periodic_box(py, cell.as_ref())?,
            })
        }
    }
}

pub(crate) fn typed_request(value: &Bound<'_, PyAny>) -> Option<PySpatialOperation> {
    if let Ok(value) = value.extract::<PyRef<'_, PyNeighborPairs>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyAtomsWithin>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    None
}

pub(crate) fn serialized_request(
    py: Python<'_>,
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PySpatialOperation>> {
    let Some(value) = config.get_item("operation")? else {
        return Ok(None);
    };
    match value.extract::<&str>()? {
        "neighbor_pairs" => Ok(Some(PyNeighborPairs::from_dict(py, config)?.operation)),
        "atoms_within" => Ok(Some(PyAtomsWithin::from_dict(py, config)?.operation)),
        _ => Ok(None),
    }
}

pub(super) fn explain_operation<'py>(
    py: Python<'py>,
    operation: &PySpatialOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation_tag(operation))?;
    result.set_item("execution", "native")?;
    result.set_item("gil_released", true)?;
    result.set_item("requires_c_contiguous", true)?;
    result.set_item("backend", backend_label(operation))?;
    result.set_item(
        "selection_contract",
        "sorted unique C-contiguous uint32 arrays",
    )?;
    result.set_item("selection_validation", "O(k) per array")?;
    result.set_item(
        "selection_materialization",
        "one owned native AtomSelection per supplied array at execution",
    )?;
    result.set_item(
        "complexity",
        "planner-dependent: O(L×R) brute force or O(N+K) indexed",
    )?;
    result.set_item(
        "materializes",
        match operation {
            PySpatialOperation::NeighborPairs { .. } => "NeighborTable columns",
            PySpatialOperation::AtomsWithin { .. } => "immutable atom-index array",
        },
    )?;
    Ok(result)
}

pub(super) fn operation_to_dict<'py>(
    py: Python<'py>,
    operation: &PySpatialOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    match operation {
        PySpatialOperation::NeighborPairs {
            positions,
            left,
            right,
            cutoff,
            options,
            cell,
        } => {
            result.set_item("operation", "neighbor_pairs")?;
            result.set_item("positions", positions.bind(py))?;
            set_optional_object(&result, "left", py, left.as_ref())?;
            set_optional_object(&result, "right", py, right.as_ref())?;
            result.set_item("cutoff", cutoff)?;
            result.set_item("options", Py::new(py, *options)?)?;
            set_optional_object(&result, "cell", py, cell.as_ref())?;
        }
        PySpatialOperation::AtomsWithin {
            positions,
            targets,
            query,
            cutoff,
            options,
            cell,
        } => {
            result.set_item("operation", "atoms_within")?;
            result.set_item("positions", positions.bind(py))?;
            result.set_item("targets", targets.bind(py))?;
            set_optional_object(&result, "query", py, query.as_ref())?;
            result.set_item("cutoff", cutoff)?;
            result.set_item("options", Py::new(py, *options)?)?;
            set_optional_object(&result, "cell", py, cell.as_ref())?;
        }
    }
    Ok(result)
}

fn operation_tag(operation: &PySpatialOperation) -> &'static str {
    match operation {
        PySpatialOperation::NeighborPairs { .. } => "neighbor_pairs",
        PySpatialOperation::AtomsWithin { .. } => "atoms_within",
    }
}

fn backend_label(operation: &PySpatialOperation) -> String {
    let options = match operation {
        PySpatialOperation::NeighborPairs { options, .. }
        | PySpatialOperation::AtomsWithin { options, .. } => options.inner(),
    };
    format!("{:?}", options.backend).to_ascii_lowercase()
}

pub(crate) fn neighbor_table_from_pairs(
    py: Python<'_>,
    pairs: Vec<molframe::NeighborPair>,
) -> PyResult<PyNeighborTable> {
    py.detach(move || PyNeighborTable::from_pairs(pairs).map_err(SpatialBindingError::Allocation))
        .map_err(SpatialBindingError::into_pyerr)
}

pub(crate) fn atom_indices_to_python(
    py: Python<'_>,
    selection: molframe::AtomSelection,
) -> Bound<'_, PyArray1<u32>> {
    let indices = py.detach(move || selection.into_iter().collect::<Vec<_>>());
    let result = numpy::ndarray::Array1::from_vec(indices).into_pyarray(py);
    let _readonly = result.readwrite().make_nonwriteable();
    result
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNeighborPairs>()?;
    module.add_class::<PyAtomsWithin>()?;
    Ok(())
}

#[cfg(test)]
#[path = "spatial_tests.rs"]
mod tests;
