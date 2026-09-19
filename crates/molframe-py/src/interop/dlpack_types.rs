//! Safe Python descriptors for the public `DLPack` ABI types.

use crate::extensions::PyExportCost;
use crate::interop::dlpack::tensor_capsule;
use crate::structure::PyStructure;
use pyo3::exceptions::{PyBufferError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use std::ptr::NonNull;

#[pyclass(name = "DLDevice", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyDLDevice {
    inner: molframe::DLDevice,
}

#[pymethods]
impl PyDLDevice {
    #[new]
    fn new(device_type: i32, device_id: i32) -> Self {
        Self {
            inner: molframe::DLDevice {
                device_type,
                device_id,
            },
        }
    }

    #[getter]
    fn device_type(slf: PyRef<'_, Self>) -> i32 {
        slf.inner.device_type
    }

    #[getter]
    fn device_id(slf: PyRef<'_, Self>) -> i32 {
        slf.inner.device_id
    }
}

#[pyclass(name = "DLDataType", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyDLDataType {
    inner: molframe::DLDataType,
}

#[pymethods]
impl PyDLDataType {
    #[new]
    fn new(code: u8, bits: u8, lanes: u16) -> Self {
        Self {
            inner: molframe::DLDataType { code, bits, lanes },
        }
    }

    #[getter]
    fn code(slf: PyRef<'_, Self>) -> u8 {
        slf.inner.code
    }

    #[getter]
    fn bits(slf: PyRef<'_, Self>) -> u8 {
        slf.inner.bits
    }

    #[getter]
    fn lanes(slf: PyRef<'_, Self>) -> u16 {
        slf.inner.lanes
    }
}

#[pyclass(name = "DLTensor", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDLTensor {
    data_address: usize,
    device: PyDLDevice,
    ndim: i32,
    dtype: PyDLDataType,
    shape: Vec<i64>,
    strides: Option<Vec<i64>>,
    byte_offset: u64,
}

#[pymethods]
impl PyDLTensor {
    #[new]
    #[pyo3(signature = (data_address, device, ndim, dtype, shape, *, strides=None, byte_offset=0))]
    fn new(
        data_address: usize,
        device: PyDLDevice,
        ndim: i32,
        dtype: PyDLDataType,
        shape: Vec<i64>,
        strides: Option<Vec<i64>>,
        byte_offset: u64,
    ) -> PyResult<Self> {
        validate_tensor_shape(ndim, shape.len(), strides.as_ref())?;
        Ok(Self {
            data_address,
            device,
            ndim,
            dtype,
            shape,
            strides,
            byte_offset,
        })
    }

    #[getter]
    const fn data_address(&self) -> usize {
        self.data_address
    }

    #[getter]
    const fn device(&self) -> PyDLDevice {
        self.device
    }

    #[getter]
    const fn ndim(&self) -> i32 {
        self.ndim
    }

    #[getter]
    const fn dtype(&self) -> PyDLDataType {
        self.dtype
    }

    #[getter]
    fn shape(&self) -> Vec<i64> {
        self.shape.clone()
    }

    #[getter]
    fn strides(&self) -> Option<Vec<i64>> {
        self.strides.clone()
    }

    #[getter]
    const fn byte_offset(&self) -> u64 {
        self.byte_offset
    }
}

#[pyclass(name = "DLManagedTensor", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDLManagedTensor {
    tensor: PyDLTensor,
    manager_context_address: usize,
    has_deleter: bool,
}

#[pymethods]
impl PyDLManagedTensor {
    #[getter]
    fn tensor(&self) -> PyDLTensor {
        self.tensor.clone()
    }

    #[getter]
    const fn manager_context_address(&self) -> usize {
        self.manager_context_address
    }

    #[getter]
    const fn has_deleter(&self) -> bool {
        self.has_deleter
    }
}

#[pyclass(name = "DlpackTensor", skip_from_py_object)]
pub(crate) struct PyDlpackTensor {
    inner: std::sync::Mutex<Option<molframe::DlpackTensor>>,
}

#[pymethods]
impl PyDlpackTensor {
    #[staticmethod]
    fn coordinates(structure: &PyStructure) -> PyResult<Self> {
        molframe::DlpackTensor::coordinates(structure.structure())
            .map(|inner| Self {
                inner: std::sync::Mutex::new(Some(inner)),
            })
            .map_err(|error| crate::errors::dlpack_error(&error))
    }

    #[getter]
    fn cost(&self) -> PyResult<PyExportCost> {
        Ok(self
            .lock()?
            .as_ref()
            .map_or(PyExportCost::Copy, |tensor| match tensor.cost() {
                molframe::ExportCost::ZeroCopy => PyExportCost::ZeroCopy,
                molframe::ExportCost::Decode => PyExportCost::Decode,
                molframe::ExportCost::Copy => PyExportCost::Copy,
            }))
    }

    #[getter]
    fn managed(&self) -> PyResult<Option<PyDLManagedTensor>> {
        Ok(self
            .lock()?
            .as_ref()
            .and_then(molframe::DlpackTensor::as_managed)
            .and_then(PyDLManagedTensor::from_native))
    }

    fn capsule<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyCapsule>> {
        let tensor = self
            .lock()?
            .take()
            .ok_or_else(|| PyBufferError::new_err("DLPack tensor has already been consumed"))?;
        tensor_capsule(py, tensor)
    }

    #[getter]
    fn consumed(&self) -> PyResult<bool> {
        Ok(self.lock()?.is_none())
    }
}

impl PyDlpackTensor {
    fn lock(&self) -> PyResult<std::sync::MutexGuard<'_, Option<molframe::DlpackTensor>>> {
        self.inner
            .lock()
            .map_err(|_| PyBufferError::new_err("DLPack tensor state is poisoned"))
    }
}

impl PyDLTensor {
    fn from_native(value: &molframe::DLTensor) -> Option<Self> {
        let ndim = value.ndim;
        if !(0..=16).contains(&ndim) {
            return None;
        }
        let count = usize::try_from(ndim).ok()?;
        let shape = if count == 0 {
            Vec::new()
        } else {
            let pointer = NonNull::new(value.shape)?;
            // SAFETY: the native DLPack producer owns `shape` for the full
            // borrow of `value`, and `ndim` was checked before constructing the
            // exact read-only slice.
            unsafe { std::slice::from_raw_parts(pointer.as_ptr(), count) }.to_vec()
        };
        let strides = NonNull::new(value.strides).map(|pointer| {
            // SAFETY: DLPack uses the same `ndim` extent for an optional strides
            // array, and the managed tensor remains alive during this copy.
            unsafe { std::slice::from_raw_parts(pointer.as_ptr(), count) }.to_vec()
        });
        Some(Self {
            data_address: value.data as usize,
            device: PyDLDevice {
                inner: value.device,
            },
            ndim,
            dtype: PyDLDataType { inner: value.dtype },
            shape,
            strides,
            byte_offset: value.byte_offset,
        })
    }
}

impl PyDLManagedTensor {
    fn from_native(value: &molframe::DLManagedTensor) -> Option<Self> {
        Some(Self {
            tensor: PyDLTensor::from_native(&value.dl_tensor)?,
            manager_context_address: value.manager_ctx as usize,
            has_deleter: value.deleter.is_some(),
        })
    }
}

fn validate_tensor_shape(ndim: i32, shape_len: usize, strides: Option<&Vec<i64>>) -> PyResult<()> {
    let ndim = usize::try_from(ndim)
        .map_err(|_| PyValueError::new_err("DLPack ndim must be non-negative"))?;
    if ndim != shape_len || strides.is_some_and(|value| value.len() != ndim) {
        return Err(PyValueError::new_err(
            "DLPack shape and strides must match ndim",
        ));
    }
    Ok(())
}

pub(crate) fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDLDevice>()?;
    module.add_class::<PyDLDataType>()?;
    module.add_class::<PyDLTensor>()?;
    module.add_class::<PyDLManagedTensor>()?;
    module.add_class::<PyDlpackTensor>()?;
    Ok(())
}
