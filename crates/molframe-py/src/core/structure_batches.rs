//! Thin Python ownership for bounded native structure batches.

use crate::core::execution::{PyBackpressure, PyBatchDemand, PyExecutionContext};
use crate::io::PyReadOptions;
use crate::structure::PyStructure;
use molframe::core::{
    Backpressure, Batch, BatchDemand, BatchLease, BatchSource, Presence, SymbolId,
};
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

type NativeLease = BatchLease<molframe::StructureBatch>;

#[pyclass(name = "StructureBatch", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureBatch(Arc<NativeLease>);

#[pymethods]
impl PyStructureBatch {
    #[getter]
    fn dataset_id(&self) -> u64 {
        self.batch().descriptor().dataset().get()
    }

    #[getter]
    fn chunk_id(&self) -> u64 {
        self.batch().descriptor().chunk().get()
    }

    #[getter]
    fn logical_start(&self) -> u64 {
        self.batch().descriptor().logical_start().get()
    }

    #[getter]
    fn rows(&self) -> usize {
        self.batch().rows()
    }

    #[getter]
    fn retained_bytes(&self) -> usize {
        self.batch().retained_bytes()
    }

    #[getter]
    fn continuity_before(&self) -> &'static str {
        continuity_name(self.batch().continuity().before)
    }

    #[getter]
    fn continuity_after(&self) -> &'static str {
        continuity_name(self.batch().continuity().after)
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat = self
            .batch()
            .positions()
            .iter()
            .flat_map(|position| position.iter().copied())
            .collect();
        let array = Array2::from_shape_vec((self.rows(), 3), flat)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?
            .into_pyarray(py);
        array.readwrite().make_nonwriteable();
        Ok(array)
    }

    #[getter]
    fn models<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        readonly(py, self.batch().models().to_vec())
    }

    #[getter]
    fn chains<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        readonly(
            py,
            self.batch()
                .chains()
                .iter()
                .map(|value| value.get())
                .collect(),
        )
    }

    #[getter]
    fn components<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        readonly(
            py,
            self.batch()
                .components()
                .iter()
                .map(|value| value.get())
                .collect(),
        )
    }

    #[getter]
    fn atoms<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        readonly(
            py,
            self.batch()
                .atoms()
                .iter()
                .map(|value| value.get())
                .collect(),
        )
    }

    #[getter]
    fn sequences<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        readonly(py, self.batch().sequences().to_vec())
    }

    #[getter]
    fn elements<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        readonly(py, self.batch().elements().to_vec())
    }

    #[getter]
    fn coordinate_presence<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        readonly(
            py,
            self.batch()
                .coordinate_presence()
                .iter()
                .copied()
                .map(presence_code)
                .collect(),
        )
    }

    #[getter]
    fn atom_site_ids<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        readonly(py, self.batch().atom_site_ids().to_vec())
    }

    #[getter]
    fn heterogens<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        readonly(py, self.batch().heterogens().to_vec())
    }

    #[getter]
    fn diagnostics(&self) -> Vec<String> {
        self.batch()
            .diagnostics()
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    fn resolve(&self, symbol: u32) -> Option<String> {
        self.batch()
            .dictionary()
            .resolve(SymbolId::from_raw(symbol))
            .map(str::to_owned)
    }
}

impl PyStructureBatch {
    fn batch(&self) -> &molframe::StructureBatch {
        self.0.batch()
    }
}

#[pyclass(name = "StructureBatchReader", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyStructureBatchReader {
    source: Mutex<molframe::StructureBatchReader>,
    context: molframe::core::ExecutionContext,
}

#[pymethods]
impl PyStructureBatchReader {
    #[pyo3(signature = (demand=None))]
    fn pull(
        &self,
        py: Python<'_>,
        demand: Option<PyBatchDemand>,
    ) -> PyResult<(PyBackpressure, Option<PyStructureBatch>)> {
        let demand = demand.map_or_else(default_demand, PyBatchDemand::native);
        let result = py.detach(|| {
            self.source
                .lock()
                .map_err(|_| PyRuntimeError::new_err("structure batch reader is poisoned"))?
                .next_batch(demand, &self.context)
                .map_err(runtime_error)
        })?;
        Ok(match result {
            Backpressure::Ready(lease) => (
                PyBackpressure::Ready,
                Some(PyStructureBatch(Arc::new(lease))),
            ),
            Backpressure::Pending => (PyBackpressure::Pending, None),
            Backpressure::Finished => (PyBackpressure::Finished, None),
        })
    }

    #[getter]
    fn context(&self) -> PyExecutionContext {
        PyExecutionContext::from_native(self.context.clone())
    }
}

#[pyfunction]
#[pyo3(signature = (path, options=None, context=None))]
pub(crate) fn open_structure_batches(
    py: Python<'_>,
    path: PathBuf,
    options: Option<&PyReadOptions>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyStructureBatchReader> {
    let options = options.map_or_else(molframe::ReadOptions::new, |value| value.0.clone());
    let context = context.map_or_else(molframe::core::ExecutionContext::default, |value| {
        value.native()
    });
    let source = py
        .detach(|| molframe::open_structure_batches(path, &options, &context))
        .map_err(runtime_error)?;
    Ok(PyStructureBatchReader {
        source: Mutex::new(source),
        context,
    })
}

#[pyfunction]
pub(crate) fn collect_structure(
    py: Python<'_>,
    source: &PyStructureBatchReader,
) -> PyResult<PyStructure> {
    let context = source.context.clone();
    let (structure, _diagnostics) = py.detach(|| {
        let mut guard = source
            .source
            .lock()
            .map_err(|_| PyRuntimeError::new_err("structure batch reader is poisoned"))?;
        molframe::collect_structure(&mut *guard, &context).map_err(runtime_error)
    })?;
    Ok(PyStructure::new(structure))
}

fn default_demand() -> BatchDemand {
    BatchDemand::new(65_536, 16 * 1024 * 1024)
}

fn continuity_name(value: molframe::ContinuityLevel) -> &'static str {
    match value {
        molframe::ContinuityLevel::None => "none",
        molframe::ContinuityLevel::Model => "model",
        molframe::ContinuityLevel::Chain => "chain",
        molframe::ContinuityLevel::Residue => "residue",
    }
}

fn presence_code(value: Presence) -> u8 {
    match value {
        Presence::Present => 0,
        Presence::Unknown => 1,
        Presence::Inapplicable => 2,
    }
}

fn readonly<T: numpy::Element>(py: Python<'_>, values: Vec<T>) -> Bound<'_, PyArray1<T>> {
    let array = values.into_pyarray(py);
    array.readwrite().make_nonwriteable();
    array
}

fn runtime_error(error: impl std::fmt::Display) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStructureBatch>()?;
    module.add_class::<PyStructureBatchReader>()?;
    module.add_function(wrap_pyfunction!(open_structure_batches, module)?)?;
    module.add_function(wrap_pyfunction!(collect_structure, module)?)
}
