//! Mechanical Python projections for native bounded MRC brick providers.

use std::path::PathBuf;

use numpy::PyArray1;
use numpy::PyArrayMethods;
use numpy::ndarray::ArrayView1;
use pyo3::prelude::*;

use crate::index::{PyChunkId, PyDatasetId};

#[pyclass(name = "MapBrickId", frozen, eq, ord, hash, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct PyMapBrickId(pub(crate) molframe::xtal::MapBrickId);

#[pymethods]
impl PyMapBrickId {
    #[new]
    fn new(value: u64) -> Self {
        Self(molframe::xtal::MapBrickId::new(value))
    }

    #[getter]
    fn value(&self) -> u64 {
        self.0.get()
    }

    fn __int__(&self) -> u64 {
        self.0.get()
    }

    fn __index__(&self) -> u64 {
        self.0.get()
    }
}

#[pyclass(name = "MrcBrickBudget", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMrcBrickBudget(pub(crate) molframe::xtal::MrcBrickBudget);

#[pymethods]
impl PyMrcBrickBudget {
    #[new]
    #[pyo3(signature = (
        *,
        max_payload_bytes=molframe::xtal::DEFAULT_MRC_BRICK_PAYLOAD_BYTES,
        max_working_set_bytes=molframe::xtal::DEFAULT_MRC_BRICK_WORKING_SET_BYTES
    ))]
    fn new(max_payload_bytes: usize, max_working_set_bytes: usize) -> Self {
        Self(molframe::xtal::MrcBrickBudget {
            max_payload_bytes,
            max_working_set_bytes,
        })
    }

    #[getter]
    fn max_payload_bytes(&self) -> usize {
        self.0.max_payload_bytes
    }

    #[getter]
    fn max_working_set_bytes(&self) -> usize {
        self.0.max_working_set_bytes
    }
}

#[pyclass(name = "MrcBrickOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMrcBrickOptions(pub(crate) molframe::xtal::MrcBrickOptions);

#[pymethods]
impl PyMrcBrickOptions {
    #[new]
    #[pyo3(signature = (*, interior=[64, 64, 64], halo=1, generation=0, budget=None))]
    fn new(
        interior: [u16; 3],
        halo: u16,
        generation: u64,
        budget: Option<PyMrcBrickBudget>,
    ) -> Self {
        Self(molframe::xtal::MrcBrickOptions {
            interior,
            halo,
            generation,
            budget: match budget {
                Some(value) => value.0,
                None => molframe::xtal::MrcBrickBudget::default(),
            },
        })
    }

    #[getter]
    fn interior(&self) -> [u16; 3] {
        self.0.interior
    }

    #[getter]
    fn halo(&self) -> u16 {
        self.0.halo
    }

    #[getter]
    fn generation(&self) -> u64 {
        self.0.generation
    }

    #[getter]
    fn budget(&self) -> PyMrcBrickBudget {
        PyMrcBrickBudget(self.0.budget)
    }
}

#[pyclass(name = "MapBrickAddress", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMapBrickAddress(molframe::xtal::MapBrickAddress);

#[pymethods]
impl PyMapBrickAddress {
    #[getter]
    fn origin(&self) -> [u64; 3] {
        self.0.origin
    }

    #[getter]
    fn mip(&self) -> u16 {
        self.0.mip
    }
}

#[pyclass(name = "MapBrickShape", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMapBrickShape(molframe::xtal::MapBrickShape);

#[pymethods]
impl PyMapBrickShape {
    #[getter]
    fn stored(&self) -> [u16; 3] {
        self.0.stored
    }

    #[getter]
    fn interior(&self) -> [u16; 3] {
        self.0.interior
    }

    #[getter]
    fn halo(&self) -> u16 {
        self.0.halo
    }

    #[getter]
    fn voxel_count(&self) -> u32 {
        self.0.voxel_count
    }
}

#[pyclass(name = "ScalarBrickMetadata", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyScalarBrickMetadata(molframe::xtal::ScalarBrickMetadata);

#[pymethods]
impl PyScalarBrickMetadata {
    #[getter]
    fn id(&self) -> PyMapBrickId {
        PyMapBrickId(self.0.id)
    }

    #[getter]
    fn address(&self) -> PyMapBrickAddress {
        PyMapBrickAddress(self.0.address)
    }

    #[getter]
    fn shape(&self) -> PyMapBrickShape {
        PyMapBrickShape(self.0.shape)
    }

    #[getter]
    fn minimum(&self) -> f32 {
        self.0.min
    }

    #[getter]
    fn maximum(&self) -> f32 {
        self.0.max
    }

    #[getter]
    fn generation(&self) -> u64 {
        self.0.generation
    }
}

#[pyclass(name = "MrcBrickDescriptor", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMrcBrickDescriptor(molframe::xtal::MrcBrickDescriptor);

#[pymethods]
impl PyMrcBrickDescriptor {
    #[getter]
    fn dataset(&self) -> PyDatasetId {
        PyDatasetId(self.0.dataset)
    }

    #[getter]
    fn chunk(&self) -> PyChunkId {
        PyChunkId(self.0.chunk)
    }

    #[getter]
    fn metadata(&self) -> PyScalarBrickMetadata {
        PyScalarBrickMetadata(self.0.metadata)
    }
}

#[pyclass(name = "ScalarBrickPayload", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyScalarBrickPayload(pub(crate) molframe::xtal::ScalarBrickPayload);

#[pymethods]
impl PyScalarBrickPayload {
    #[getter]
    fn descriptor(&self) -> PyMrcBrickDescriptor {
        PyMrcBrickDescriptor(self.0.descriptor())
    }

    #[getter]
    fn actual_range(&self) -> [f32; 2] {
        self.0.actual_range()
    }

    #[getter]
    fn payload_bytes(&self) -> usize {
        self.0.payload_bytes()
    }

    #[getter]
    fn values<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let borrowed = slf.borrow();
        let values = borrowed.0.values();
        // SAFETY: the view points into the immutable contiguous `Vec<f32>`
        // retained by `slf`, which becomes NumPy's base owner below.
        let view = unsafe { ArrayView1::from_shape_ptr(values.len(), values.as_ptr()) };
        // SAFETY: NumPy retains `slf` for at least the complete view lifetime.
        let array = unsafe { PyArray1::borrow_from_array(&view, slf.clone().into_any()) };
        let _readonly = array.readwrite().make_nonwriteable();
        array
    }
}

#[pyclass(name = "MrcBrickProvider", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyMrcBrickProvider(molframe::xtal::MrcBrickProvider<std::fs::File>);

#[pymethods]
impl PyMrcBrickProvider {
    #[new]
    #[pyo3(signature = (path, dataset, first_chunk, first_brick, *, options=None))]
    fn new(
        path: PathBuf,
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        first_brick: PyMapBrickId,
        options: Option<PyMrcBrickOptions>,
    ) -> PyResult<Self> {
        let options = match options {
            Some(value) => value.0,
            None => molframe::xtal::MrcBrickOptions::default(),
        };
        molframe::xtal::MrcBrickProvider::open(
            path,
            dataset.0,
            first_chunk.0,
            first_brick.0,
            options,
        )
        .map(Self)
        .map_err(brick_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetId {
        PyDatasetId(self.0.dataset())
    }

    #[getter]
    fn logical_extent(&self) -> [u64; 3] {
        self.0.logical_extent()
    }

    #[getter]
    fn brick_count(&self) -> u64 {
        self.0.brick_count()
    }

    fn descriptor(&self, id: PyMapBrickId) -> PyResult<PyMrcBrickDescriptor> {
        self.0
            .descriptor(id.0)
            .map(PyMrcBrickDescriptor)
            .map_err(brick_error)
    }

    fn read_brick(&mut self, id: PyMapBrickId) -> PyResult<PyScalarBrickPayload> {
        self.0
            .read_brick(id.0)
            .map(PyScalarBrickPayload)
            .map_err(brick_error)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_MRC_BRICK_PAYLOAD_BYTES",
        molframe::xtal::DEFAULT_MRC_BRICK_PAYLOAD_BYTES,
    )?;
    module.add(
        "DEFAULT_MRC_BRICK_WORKING_SET_BYTES",
        molframe::xtal::DEFAULT_MRC_BRICK_WORKING_SET_BYTES,
    )?;
    module.add_class::<PyMapBrickId>()?;
    module.add_class::<PyMrcBrickBudget>()?;
    module.add_class::<PyMrcBrickOptions>()?;
    module.add_class::<PyMapBrickAddress>()?;
    module.add_class::<PyMapBrickShape>()?;
    module.add_class::<PyScalarBrickMetadata>()?;
    module.add_class::<PyMrcBrickDescriptor>()?;
    module.add_class::<PyScalarBrickPayload>()?;
    module.add_class::<PyMrcBrickProvider>()
}

fn brick_error(error: impl std::fmt::Display) -> PyErr {
    crate::errors::MrcBrickError::new_err(error.to_string())
}

#[cfg(test)]
#[path = "bricks_tests.rs"]
mod tests;
