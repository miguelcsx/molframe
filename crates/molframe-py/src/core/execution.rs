//! Declarative Python controls for bounded native execution.

use molframe::core::SourceBytes;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::path::PathBuf;
use std::sync::Mutex;

#[pyclass(name = "MemoryBudget", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMemoryBudget(molframe::core::MemoryBudget);

#[pymethods]
impl PyMemoryBudget {
    #[new]
    #[pyo3(signature = (bytes=molframe::core::execution::DEFAULT_MEMORY_BUDGET_BYTES))]
    fn new(bytes: usize) -> PyResult<Self> {
        molframe::core::MemoryBudget::new(bytes)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn bytes(&self) -> usize {
        self.0.bytes()
    }
}

#[pyclass(name = "CancellationToken", frozen, from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyCancellationToken(molframe::core::CancellationToken);

#[pymethods]
impl PyCancellationToken {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn cancel(&self) {
        self.0.cancel();
    }

    #[getter]
    fn is_cancelled(&self) -> bool {
        self.0.is_cancelled()
    }

    fn check(&self) -> PyResult<()> {
        self.0.check().map_err(runtime_error)
    }
}

#[pyclass(name = "ScratchPolicy", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyScratchPolicy(molframe::core::ScratchPolicy);

#[pymethods]
impl PyScratchPolicy {
    #[new]
    #[pyo3(signature = (max_bytes=8_000_000))]
    fn new(max_bytes: usize) -> Self {
        Self(molframe::core::ScratchPolicy::new(max_bytes))
    }

    #[getter]
    fn max_bytes(&self) -> usize {
        self.0.max_bytes()
    }
}

#[pyclass(name = "TempStoragePolicy", frozen, from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyTempStoragePolicy(molframe::core::TempStoragePolicy);

#[pymethods]
impl PyTempStoragePolicy {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    #[staticmethod]
    fn disabled() -> Self {
        Self::default()
    }

    #[staticmethod]
    fn directory(root: PathBuf, max_bytes: u64) -> Self {
        Self(molframe::core::TempStoragePolicy::directory(
            root, max_bytes,
        ))
    }

    #[getter]
    fn root(&self) -> Option<PathBuf> {
        self.0.root().map(std::path::Path::to_path_buf)
    }

    #[getter]
    fn max_bytes(&self) -> u64 {
        self.0.max_bytes()
    }
}

#[pyclass(name = "ExecutionContext", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyExecutionContext(molframe::core::ExecutionContext);

impl PyExecutionContext {
    pub(crate) fn native(&self) -> molframe::core::ExecutionContext {
        self.0.clone()
    }

    pub(crate) fn from_native(value: molframe::core::ExecutionContext) -> Self {
        Self(value)
    }
}

pub(crate) fn default_context() -> molframe::core::ExecutionContext {
    molframe::core::ExecutionContext::default()
}

#[pyclass(name = "WindowedFile", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyWindowedFile(Mutex<molframe::core::WindowedFile>);

#[pymethods]
impl PyWindowedFile {
    #[new]
    fn new(path: PathBuf, max_window_bytes: usize, context: &PyExecutionContext) -> PyResult<Self> {
        molframe::core::WindowedFile::open(path, max_window_bytes, &context.0)
            .map(|source| Self(Mutex::new(source)))
            .map_err(runtime_error)
    }

    #[getter]
    fn capacity(&self) -> PyResult<usize> {
        self.0
            .lock()
            .map(|source| source.capacity())
            .map_err(|_| PyRuntimeError::new_err("windowed file state is poisoned"))
    }

    #[getter]
    fn length(&self) -> PyResult<u64> {
        self.0
            .lock()
            .map_err(|_| PyRuntimeError::new_err("windowed file state is poisoned"))?
            .len_hint()
            .ok_or_else(|| PyRuntimeError::new_err("windowed file length is unavailable"))
    }

    fn window<'py>(
        &self,
        py: Python<'py>,
        start: u64,
        length: usize,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let bytes = py.detach(|| {
            let mut source = self
                .0
                .lock()
                .map_err(|_| PyRuntimeError::new_err("windowed file state is poisoned"))?;
            source
                .window(start, length)
                .map(|window| window.bytes().to_vec())
                .map_err(runtime_error)
        })?;
        Ok(PyBytes::new(py, &bytes))
    }
}

#[pymethods]
impl PyExecutionContext {
    #[new]
    #[pyo3(signature = (*, memory_budget=None, cancellation=None, scratch_policy=None, temp_storage_policy=None, worker_budget=None))]
    fn new(
        memory_budget: Option<PyMemoryBudget>,
        cancellation: Option<PyCancellationToken>,
        scratch_policy: Option<PyScratchPolicy>,
        temp_storage_policy: Option<PyTempStoragePolicy>,
        worker_budget: Option<usize>,
    ) -> PyResult<Self> {
        let mut builder = molframe::core::ExecutionContext::builder();
        if let Some(value) = memory_budget {
            builder = builder.memory_budget(value.0);
        }
        if let Some(value) = cancellation {
            builder = builder.cancellation(value.0);
        }
        if let Some(value) = scratch_policy {
            builder = builder.scratch_policy(value.0);
        }
        if let Some(value) = temp_storage_policy {
            builder = builder.temp_storage_policy(value.0);
        }
        if let Some(value) = worker_budget {
            builder = builder.worker_budget(value);
        }
        builder.build().map(Self).map_err(value_error)
    }

    #[getter]
    fn memory_budget(&self) -> PyMemoryBudget {
        PyMemoryBudget(self.0.memory_budget())
    }

    #[getter]
    fn worker_budget(&self) -> usize {
        self.0.worker_budget()
    }

    #[getter]
    fn cancellation(&self) -> PyCancellationToken {
        PyCancellationToken(self.0.cancellation().clone())
    }

    #[getter]
    fn scratch_policy(&self) -> PyScratchPolicy {
        PyScratchPolicy(self.0.scratch_policy())
    }

    #[getter]
    fn temp_storage_policy(&self) -> PyTempStoragePolicy {
        PyTempStoragePolicy(self.0.temp_storage_policy().clone())
    }

    fn reserve(&self, bytes: usize) -> PyResult<PyMemoryLease> {
        self.0
            .try_reserve(bytes)
            .map(PyMemoryLease)
            .map_err(runtime_error)
    }

    fn create_spill_file(&self, key: &str) -> PyResult<PySpillFile> {
        self.0
            .create_spill_file(key)
            .map(|file| PySpillFile(Some(file)))
            .map_err(runtime_error)
    }

    #[getter]
    fn spill_bytes(&self) -> u64 {
        self.0.spill_bytes()
    }

    #[getter]
    fn total_spilled_bytes(&self) -> u64 {
        self.0.total_spilled_bytes()
    }

    #[getter]
    fn reserved_bytes(&self) -> usize {
        self.0.reserved_bytes()
    }

    #[getter]
    fn peak_reserved_bytes(&self) -> usize {
        self.0.peak_reserved_bytes()
    }

    #[getter]
    fn live_batches(&self) -> usize {
        self.0.live_batches()
    }

    #[getter]
    fn peak_live_batches(&self) -> usize {
        self.0.peak_live_batches()
    }
}

#[pyclass(name = "MemoryLease", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyMemoryLease(molframe::core::MemoryReservation);

#[pymethods]
impl PyMemoryLease {
    #[getter]
    fn bytes(&self) -> usize {
        self.0.bytes()
    }
}

#[pyclass(name = "SpillFile", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PySpillFile(Option<molframe::core::SpillFile>);

#[pymethods]
impl PySpillFile {
    fn append(&mut self, payload: &Bound<'_, PyBytes>) -> PyResult<()> {
        let Some(file) = &mut self.0 else {
            return Err(PyRuntimeError::new_err("spill file is already finished"));
        };
        file.append(payload.as_bytes()).map_err(runtime_error)
    }

    fn finish(&mut self) -> PyResult<PySpillArtifact> {
        let Some(file) = self.0.take() else {
            return Err(PyRuntimeError::new_err("spill file is already finished"));
        };
        file.finish().map(PySpillArtifact).map_err(runtime_error)
    }
}

#[pyclass(name = "SpillArtifact", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PySpillArtifact(molframe::core::SpillArtifact);

#[pymethods]
impl PySpillArtifact {
    #[getter]
    fn bytes(&self) -> u64 {
        self.0.bytes()
    }

    #[getter]
    fn records(&self) -> u64 {
        self.0.records()
    }

    #[getter]
    fn path(&self) -> PathBuf {
        self.0.path().to_owned()
    }

    fn reader(&self) -> PyResult<PySpillReader> {
        self.0
            .reader()
            .map(|reader| PySpillReader {
                reader,
                buffer: Vec::new(),
            })
            .map_err(runtime_error)
    }
}

#[pyclass(name = "SpillReader", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PySpillReader {
    reader: molframe::core::SpillReader,
    buffer: Vec<u8>,
}

#[pymethods]
impl PySpillReader {
    #[pyo3(signature = (max_bytes=8_000_000))]
    fn read_next<'py>(
        &mut self,
        py: Python<'py>,
        max_bytes: usize,
    ) -> PyResult<Option<Bound<'py, PyBytes>>> {
        match self
            .reader
            .read_next(&mut self.buffer, max_bytes)
            .map_err(runtime_error)?
        {
            Some(_length) => Ok(Some(PyBytes::new(py, &self.buffer))),
            None => Ok(None),
        }
    }
}

#[pyclass(name = "BatchDemand", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBatchDemand(molframe::core::BatchDemand);

impl PyBatchDemand {
    pub(crate) const fn native(self) -> molframe::core::BatchDemand {
        self.0
    }
}

#[pymethods]
impl PyBatchDemand {
    #[new]
    fn new(max_rows: usize, max_bytes: usize) -> Self {
        Self(molframe::core::BatchDemand::new(max_rows, max_bytes))
    }

    #[getter]
    fn max_rows(&self) -> usize {
        self.0.max_rows
    }

    #[getter]
    fn max_bytes(&self) -> usize {
        self.0.max_bytes
    }

    #[getter]
    fn can_accept_work(&self) -> bool {
        self.0.can_accept_work()
    }

    fn accepts(&self, rows: usize, retained_bytes: usize) -> bool {
        rows <= self.0.max_rows && retained_bytes <= self.0.max_bytes
    }
}

#[pyclass(name = "Backpressure", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyBackpressure {
    Ready,
    Pending,
    Finished,
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn runtime_error(error: impl std::fmt::Display) -> PyErr {
    PyRuntimeError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_MEMORY_BUDGET_BYTES",
        molframe::core::execution::DEFAULT_MEMORY_BUDGET_BYTES,
    )?;
    module.add_class::<PyMemoryBudget>()?;
    module.add_class::<PyCancellationToken>()?;
    module.add_class::<PyScratchPolicy>()?;
    module.add_class::<PyTempStoragePolicy>()?;
    module.add_class::<PyExecutionContext>()?;
    module.add_class::<PyWindowedFile>()?;
    module.add_class::<PyMemoryLease>()?;
    module.add_class::<PySpillFile>()?;
    module.add_class::<PySpillArtifact>()?;
    module.add_class::<PySpillReader>()?;
    module.add_class::<PyBatchDemand>()?;
    module.add_class::<PyBackpressure>()
}
