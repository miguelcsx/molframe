//! Compact core columns and validity masks.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "BitVec", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBitVec(pub(crate) molframe::core::BitVec);

#[pymethods]
impl PyBitVec {
    #[new]
    fn new() -> Self {
        Self(molframe::core::BitVec::new())
    }

    #[staticmethod]
    fn repeat(value: bool, length: u32) -> Self {
        Self(molframe::core::BitVec::repeat(value, length))
    }

    #[staticmethod]
    fn with_capacity(length: u32) -> Self {
        Self(molframe::core::BitVec::with_capacity(length))
    }

    #[staticmethod]
    fn from_values(values: Vec<bool>) -> PyResult<Self> {
        molframe::core::BitVec::try_from_iter(values)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err("bit vector length exceeds the native limit"))
    }

    #[staticmethod]
    fn try_from_iter(values: Vec<bool>) -> PyResult<Self> {
        Self::from_values(values)
    }

    #[getter]
    fn len(&self) -> u32 {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn try_push(&mut self, value: bool) -> PyResult<()> {
        self.0
            .try_push(value)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn get(&self, position: u32) -> Option<bool> {
        self.0.get(position)
    }

    fn test(&self, position: u32) -> bool {
        self.0.test(position)
    }

    fn set(&mut self, position: u32, value: bool) {
        self.0.set(position, value);
    }

    fn count_ones(&self) -> u32 {
        self.0.count_ones()
    }

    fn all(&self) -> bool {
        self.0.all()
    }

    fn none(&self) -> bool {
        self.0.none()
    }

    fn ones(&self) -> Vec<u32> {
        self.0.ones().collect()
    }

    fn intersect_with(&mut self, other: &PyBitVec) {
        self.0.intersect_with(&other.0);
    }

    fn union_with(&mut self, other: &PyBitVec) {
        self.0.union_with(&other.0);
    }

    fn invert(&mut self) {
        self.0.invert();
    }
}

#[pyclass(name = "Presence", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPresence {
    Present,
    Unknown,
    Inapplicable,
}

impl From<PyPresence> for molframe::core::Presence {
    fn from(value: PyPresence) -> Self {
        match value {
            PyPresence::Present => Self::Present,
            PyPresence::Unknown => Self::Unknown,
            PyPresence::Inapplicable => Self::Inapplicable,
        }
    }
}

impl From<molframe::core::Presence> for PyPresence {
    fn from(value: molframe::core::Presence) -> Self {
        match value {
            molframe::core::Presence::Present => Self::Present,
            molframe::core::Presence::Unknown => Self::Unknown,
            molframe::core::Presence::Inapplicable => Self::Inapplicable,
        }
    }
}

#[pymethods]
impl PyPresence {
    fn is_present(&self) -> bool {
        *self == Self::Present
    }
}

#[pyclass(name = "ValidityMask", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyValidityMask(pub(crate) molframe::core::ValidityMask);

#[pymethods]
impl PyValidityMask {
    #[staticmethod]
    fn all_present(length: u32) -> Self {
        Self(molframe::core::ValidityMask::all_present(length))
    }

    #[staticmethod]
    fn from_values(values: Vec<PyPresence>) -> PyResult<Self> {
        let values = values.into_iter().map(Into::into).collect::<Vec<_>>();
        molframe::core::ValidityMask::try_from_iter(values)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err("validity mask length exceeds the native limit"))
    }

    #[staticmethod]
    fn try_from_iter(values: Vec<PyPresence>) -> PyResult<Self> {
        Self::from_values(values)
    }

    #[getter]
    fn len(&self) -> u32 {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn is_all_present(&self) -> bool {
        self.0.is_all_present()
    }

    fn get(&self, position: u32) -> PyPresence {
        self.0.get(position).into()
    }

    fn set(&mut self, position: u32, presence: PyPresence) {
        self.0.set(position, presence.into());
    }

    fn present_count(&self) -> u32 {
        self.0.present_count()
    }

    fn compact(&mut self) {
        self.0.compact();
    }
}

#[pyfunction]
pub(crate) fn bit_width(py: Python<'_>, maximum: u64) -> u8 {
    py.detach(move || -> u8 { molframe::core::column::bit_width(maximum) })
}

#[pyfunction(name = "pack")]
pub(crate) fn pack_bits(values: Vec<u64>, width: u8) -> PyResult<Vec<u8>> {
    molframe::core::column::pack(&values, width)
        .ok_or_else(|| PyValueError::new_err("packed bit column exceeds the native limit"))
}

#[pyfunction(name = "unpack_one")]
pub(crate) fn unpack_one_bits(data: Vec<u8>, width: u8, index: u32) -> Option<u64> {
    molframe::core::column::unpack_one(&data, width, index)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBitVec>()?;
    module.add_class::<PyPresence>()?;
    module.add_class::<PyValidityMask>()?;
    module.add_function(wrap_pyfunction!(bit_width, module)?)?;
    module.add_function(wrap_pyfunction!(pack_bits, module)?)?;
    module.add_function(wrap_pyfunction!(unpack_one_bits, module)?)?;
    Ok(())
}
