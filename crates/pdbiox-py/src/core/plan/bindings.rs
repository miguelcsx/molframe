//! Declarative native operations and deterministic Rust-side composition.

#[path = "chemistry.rs"]
mod chemistry;
#[path = "coordinates.rs"]
mod coordinates;
#[path = "dispatch.rs"]
mod dispatch;
#[path = "geometry.rs"]
mod geometry;
#[path = "lowering.rs"]
mod lowering;
#[path = "physical.rs"]
mod physical;
#[path = "results.rs"]
mod results;
#[path = "selection.rs"]
mod selection;
#[path = "spatial.rs"]
mod spatial;
#[path = "structure.rs"]
mod structure;
#[path = "surface.rs"]
mod surface;
#[path = "trajectory.rs"]
mod trajectory;

use crate::geometry::borrowed_coordinates;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use coordinates::{CoordinateMetric, PyComparison};
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyString};
use std::collections::{BTreeMap, HashMap};

use dispatch::operation_from_value;
use lowering::{lower_operation, plan_error};

pub(crate) use coordinates::{PyContacts, PyGdtHa, PyGdtTs, PyLddt, PyRmsd, PyTmScore};
enum Operation {
    Contacts(PyContacts),
    Rmsd(PyRmsd),
    Comparison(PyComparison),
    Selection(selection::PySelectQuery),
    BasePairs(structure::PyBasePairs),
    Structure(pdbiox::StructureRequest),
    BondInference(chemistry::PyInferBonds),
    Geometry(geometry::PyGeometryOperation),
    Spatial(spatial::PySpatialOperation),
    Physical(physical::PyPhysicalOperation),
    Surface(surface::PySurfaceOperation),
    Trajectory(trajectory::PyTrajectoryOperation),
}

#[pyclass(name = "Plan")]
pub(crate) struct PyPlan {
    operations: Vec<(String, Operation)>,
}

#[pyclass(name = "PlanResult", frozen)]
pub(crate) struct PyPlanResult {
    #[pyo3(get)]
    results: Py<PyDict>,
    #[pyo3(get)]
    cached_index_count: usize,
}

#[pymethods]
impl PyPlan {
    #[new]
    #[pyo3(signature = (**operations))]
    fn new(operations: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut entries = BTreeMap::new();
        let Some(operations) = operations else {
            return Err(PyValueError::new_err(
                "Plan requires at least one operation",
            ));
        };
        for (key, value) in operations.iter() {
            let name = key
                .cast::<PyString>()
                .map_err(|_| PyTypeError::new_err("Plan operation ids must be strings"))?
                .to_str()?
                .to_owned();
            if entries.contains_key(&name) {
                return Err(PyKeyError::new_err(format!(
                    "duplicate operation id: {name}"
                )));
            }
            entries.insert(name.clone(), operation_from_value(value, &name)?);
        }
        Ok(Self {
            operations: entries.into_iter().collect(),
        })
    }

    #[pyo3(signature = (structure=None, *, context=None))]
    fn execute<'py>(
        &self,
        py: Python<'py>,
        structure: Option<&'py PyStructure>,
        context: Option<&crate::core::execution::PyExecutionContext>,
    ) -> PyResult<Bound<'py, PyDict>> {
        results::plan_result_dict(py, execute_native(self, py, structure, context)?, structure)
    }

    #[pyo3(signature = (structure=None, *, context=None))]
    fn execute_report<'py>(
        &self,
        py: Python<'py>,
        structure: Option<&'py PyStructure>,
        context: Option<&crate::core::execution::PyExecutionContext>,
    ) -> PyResult<PyPlanResult> {
        let native = execute_native(self, py, structure, context)?;
        let cached_index_count = native.cached_index_count;
        let results = results::plan_result_dict(py, native, structure)?.unbind();
        Ok(PyPlanResult {
            results,
            cached_index_count,
        })
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("execution", "rust")?;
        result.set_item("deterministic_order", true)?;
        result.set_item("operation_count", self.operations.len())?;
        let operations = PyDict::new(py);
        for (name, operation) in &self.operations {
            let explanation = match operation {
                Operation::Contacts(value) => value.explain(py)?,
                Operation::Rmsd(value) => value.explain(py)?,
                Operation::Comparison(value) => coordinates::comparison_explain(py, value.metric)?,
                Operation::Selection(value) => selection::explain(py, value)?.unbind().into_any(),
                Operation::BasePairs(value) => value.explain(py)?,
                Operation::Structure(value) => structure::explain_request(py, value)?,
                Operation::BondInference(value) => value.explain_native(py)?.unbind().into_any(),
                Operation::Geometry(value) => {
                    geometry::explain_operation(py, value)?.unbind().into_any()
                }
                Operation::Spatial(value) => {
                    spatial::explain_operation(py, value)?.unbind().into_any()
                }
                Operation::Physical(value) => {
                    physical::explain_operation(py, value)?.unbind().into_any()
                }
                Operation::Surface(value) => {
                    surface::explain_operation(py, value)?.unbind().into_any()
                }
                Operation::Trajectory(value) => trajectory::explain_operation(py, value)?
                    .unbind()
                    .into_any(),
            };
            operations.set_item(name, explanation)?;
        }
        result.set_item("operations", operations)?;
        Ok(result)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        let operations = PyDict::new(py);
        for (name, operation) in &self.operations {
            let config = match operation {
                Operation::Contacts(value) => value.to_dict(py)?,
                Operation::Rmsd(value) => value.to_dict(py)?,
                Operation::Comparison(value) => coordinates::comparison_to_dict(
                    py,
                    &value.mobile,
                    &value.reference,
                    value.metric,
                )?,
                Operation::Selection(value) => selection::to_dict(py, value)?.unbind().into_any(),
                Operation::BasePairs(value) => value.to_dict(py)?,
                Operation::Structure(value) => structure::dict_request(py, value)?,
                Operation::BondInference(value) => value.to_dict_native(py)?.unbind().into_any(),
                Operation::Geometry(value) => {
                    geometry::operation_to_dict(py, value)?.unbind().into_any()
                }
                Operation::Spatial(value) => {
                    spatial::operation_to_dict(py, value)?.unbind().into_any()
                }
                Operation::Physical(value) => {
                    physical::operation_to_dict(py, value)?.unbind().into_any()
                }
                Operation::Surface(value) => {
                    surface::operation_to_dict(py, value)?.unbind().into_any()
                }
                Operation::Trajectory(value) => trajectory::operation_to_dict(py, value)?
                    .unbind()
                    .into_any(),
            };
            operations.set_item(name, config.bind(py))?;
        }
        result.set_item("operations", operations)?;
        Ok(result)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operations = required(config, "operations")?.cast_into::<PyDict>()?;
        Self::new(Some(&operations))
    }

    fn __repr__(&self) -> String {
        let names = self
            .operations
            .iter()
            .map(|(name, _)| name.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        format!("Plan({names})")
    }
}

/// Lowers Python operation objects into the native facade plan. Only the
/// operation graph is assembled while the GIL is held; scientific data-plane
/// work runs after `detach` in Rust.
fn execute_native(
    plan: &PyPlan,
    py: Python<'_>,
    structure: Option<&PyStructure>,
    context: Option<&crate::core::execution::PyExecutionContext>,
) -> PyResult<pdbiox::PlanResult> {
    let mut native = pdbiox::Plan::new();
    let mut arrays = Vec::new();
    let mut scalar_arrays = Vec::new();
    let mut float_arrays = Vec::new();
    let mut mask_arrays = Vec::new();
    let mut frame_arrays = Vec::new();
    let mut index_arrays = Vec::new();
    let mut coordinate_slots = HashMap::new();
    let mut scalar_slots = HashMap::new();
    let mut float_slots = HashMap::new();
    let mut mask_slots = HashMap::new();
    let mut frame_slots = HashMap::new();
    let mut index_slots = HashMap::new();

    for (name, operation) in &plan.operations {
        lower_operation(
            &mut native,
            name,
            operation,
            py,
            &mut arrays,
            &mut scalar_arrays,
            &mut float_arrays,
            &mut mask_arrays,
            &mut frame_arrays,
            &mut index_arrays,
            &mut coordinate_slots,
            &mut scalar_slots,
            &mut float_slots,
            &mut mask_slots,
            &mut frame_slots,
            &mut index_slots,
        )?;
    }

    execute_lowered(
        py,
        native,
        structure,
        arrays,
        scalar_arrays,
        float_arrays,
        mask_arrays,
        frame_arrays,
        index_arrays,
        context,
    )
}

fn execute_lowered<'py>(
    py: Python<'py>,
    native: pdbiox::Plan,
    structure: Option<&PyStructure>,
    arrays: Vec<PyReadonlyArray2<'py, f32>>,
    scalar_arrays: Vec<numpy::PyReadonlyArray1<'py, f64>>,
    float_arrays: Vec<numpy::PyReadonlyArray1<'py, f32>>,
    mask_arrays: Vec<numpy::PyReadonlyArray1<'py, bool>>,
    frame_arrays: Vec<numpy::PyReadonlyArray3<'py, f32>>,
    index_arrays: Vec<PyReadonlyArray1<'py, usize>>,
    context: Option<&crate::core::execution::PyExecutionContext>,
) -> PyResult<pdbiox::PlanResult> {
    let coordinates = arrays
        .iter()
        .map(borrowed_coordinates)
        .collect::<PyResult<Vec<_>>>()?;
    let inputs = coordinates
        .iter()
        .map(|positions| pdbiox::CoordinateInput { positions })
        .collect::<Vec<_>>();
    let scalar_values = scalar_arrays
        .iter()
        .map(|array| {
            array.as_slice().map_err(|_| {
                PyValueError::new_err(
                    "scalar inputs must be C-contiguous; pass copy=True explicitly to materialise them",
                )
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let scalars = scalar_values
        .iter()
        .map(|values| pdbiox::ScalarInput { values })
        .collect::<Vec<_>>();
    let float_values = float_arrays
        .iter()
        .map(|array| {
            array.as_slice().map_err(|_| {
                PyValueError::new_err(
                    "floating-point inputs must be C-contiguous; pass copy=True explicitly to materialise them",
                )
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let floats = float_values
        .iter()
        .map(|values| pdbiox::FloatInput { values })
        .collect::<Vec<_>>();
    let mask_values = mask_arrays
        .iter()
        .map(|array| {
            array.as_slice().map_err(|_| {
                PyValueError::new_err(
                    "surface mask inputs must be C-contiguous; pass copy=True explicitly to materialise them",
                )
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let masks = mask_values
        .iter()
        .map(|values| pdbiox::MaskInput { values })
        .collect::<Vec<_>>();
    let frame_values = frame_inputs(&frame_arrays)?;
    let index_values = index_arrays
        .iter()
        .map(|array| {
            array.as_slice().map_err(|_| {
                PyValueError::new_err(
                    "atom-index inputs must be C-contiguous; call numpy.ascontiguousarray explicitly",
                )
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let indices = index_values
        .iter()
        .map(|indices| pdbiox::IndexInput { indices })
        .collect::<Vec<_>>();
    let retained_structure = structure.map(|value| value.structure().clone());
    let execution = context.map_or_else(pdbiox::core::ExecutionContext::default, |value| {
        value.native()
    });
    py.detach(|| {
        native
            .execute(
                pdbiox::PlanInput {
                    structure: retained_structure.as_ref(),
                    arrays: &inputs,
                    scalars: &scalars,
                    floats: &floats,
                    masks: &masks,
                    frames: &frame_values,
                    indices: &indices,
                },
                &execution,
            )
            .map_err(plan_error)
    })
}

fn frame_inputs<'py>(
    arrays: &'py [numpy::PyReadonlyArray3<'py, f32>],
) -> PyResult<Vec<pdbiox::FrameInput<'py>>> {
    arrays
        .iter()
        .map(|array| crate::intrinsic::borrowed_frame_input(array))
        .collect()
}

fn execute_comparison(
    py: Python<'_>,
    mobile: &Py<PyAny>,
    reference: &Py<PyAny>,
    metric: CoordinateMetric,
) -> PyResult<f64> {
    let plan = PyPlan {
        operations: vec![(
            "result".to_owned(),
            Operation::Comparison(PyComparison {
                mobile: mobile.clone_ref(py),
                reference: reference.clone_ref(py),
                metric,
            }),
        )],
    };
    let result = execute_native(&plan, py, None, None)?;
    let Some(entry) = result.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native comparison plan returned no result",
        ));
    };
    match entry.value {
        pdbiox::PlanValue::Comparison(value) => Ok(value.value),
        _ => Err(PyValueError::new_err(
            "native comparison plan returned an incompatible result",
        )),
    }
}

fn required<'py>(config: &'py Bound<'py, PyDict>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    config
        .get_item(name)?
        .ok_or_else(|| PyKeyError::new_err(format!("missing operation field: {name}")))
}

fn optional<'py>(
    config: &'py Bound<'py, PyDict>,
    name: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    config.get_item(name)
}

impl PyAnalysisPolicy {
    fn new_default() -> Self {
        Self {
            inner: pdbiox::AnalysisPolicy::default(),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    coordinates::register(module)?;
    chemistry::register(module)?;
    geometry::register(module)?;
    selection::register(module)?;
    spatial::register(module)?;
    physical::register(module)?;
    surface::register(module)?;
    structure::register(module)?;
    trajectory::register(module)?;
    module.add_class::<PyPlan>()?;
    module.add_class::<PyPlanResult>()?;
    Ok(())
}
