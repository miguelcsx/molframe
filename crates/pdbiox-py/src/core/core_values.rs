//! Small, allocation-free core values exposed at the Python boundary.

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

create_exception!(_native, DictionaryFull, PyValueError);

#[pyclass(name = "Position", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyPosition(pub(crate) pdbiox::Position);

#[pymethods]
impl PyPosition {
    #[new]
    fn new(byte_offset: u32, line: u32, column: u32) -> Self {
        Self(pdbiox::Position::new(byte_offset, line, column))
    }

    #[classattr]
    #[pyo3(name = "START")]
    fn start() -> Self {
        Self(pdbiox::Position::START)
    }

    #[getter]
    fn byte_offset(&self) -> u32 {
        self.0.byte_offset
    }

    #[getter]
    fn line(&self) -> u32 {
        self.0.line
    }

    #[getter]
    fn column(&self) -> u32 {
        self.0.column
    }

    fn advance(&self, byte: u8) -> Option<Self> {
        self.0.advance(byte).map(Self)
    }

    fn __repr__(&self) -> String {
        format!(
            "Position({}, {}, {})",
            self.0.byte_offset, self.0.line, self.0.column
        )
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }
}

#[pyclass(name = "ByteSpan", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyByteSpan(pub(crate) pdbiox::ByteSpan);

#[pymethods]
impl PyByteSpan {
    #[new]
    fn new(start: PyPosition, end: u32) -> Self {
        Self(pdbiox::ByteSpan::new(start.0, end))
    }

    #[staticmethod]
    fn empty(at: PyPosition) -> Self {
        Self(pdbiox::ByteSpan::empty(at.0))
    }

    #[getter]
    fn start(&self) -> PyPosition {
        PyPosition(self.0.start)
    }

    #[getter]
    fn end(&self) -> u32 {
        self.0.end
    }

    fn len(&self) -> Option<u32> {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn slice<'py>(
        &self,
        py: Python<'py>,
        source: &Bound<'_, PyBytes>,
    ) -> Option<Bound<'py, PyBytes>> {
        self.0
            .slice(source.as_bytes())
            .map(|value| PyBytes::new(py, value))
    }

    fn __repr__(&self) -> String {
        format!("ByteSpan(start={:?}, end={})", self.start(), self.0.end)
    }
}

#[pyclass(name = "CoordinateGeneration", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct PyCoordinateGeneration(pub(crate) pdbiox::CoordinateGeneration);

#[pymethods]
impl PyCoordinateGeneration {
    #[new]
    fn new(value: u64) -> Self {
        Self(pdbiox::CoordinateGeneration::from_raw(value))
    }

    #[staticmethod]
    fn from_raw(value: u64) -> Self {
        Self(pdbiox::CoordinateGeneration::from_raw(value))
    }

    #[classattr]
    #[pyo3(name = "INITIAL")]
    fn initial() -> Self {
        Self(pdbiox::CoordinateGeneration::INITIAL)
    }

    #[getter]
    fn value(&self) -> u64 {
        self.0.get()
    }

    fn get(&self) -> u64 {
        self.0.get()
    }

    fn next(&self) -> Option<Self> {
        self.0.next().map(Self)
    }

    fn __int__(&self) -> u64 {
        self.0.get()
    }

    fn __index__(&self) -> u64 {
        self.0.get()
    }
}

#[pyclass(name = "Aabb", eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PyAabb(pub(crate) pdbiox::Aabb);

#[pymethods]
impl PyAabb {
    #[new]
    fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self(pdbiox::Aabb { min, max })
    }

    #[staticmethod]
    fn empty() -> Self {
        Self(pdbiox::Aabb::EMPTY)
    }

    #[classattr]
    #[pyo3(name = "EMPTY")]
    fn empty_value() -> Self {
        Self(pdbiox::Aabb::EMPTY)
    }

    #[getter]
    fn min(&self) -> [f32; 3] {
        self.0.min
    }

    #[getter]
    fn max(&self) -> [f32; 3] {
        self.0.max
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn extend(&mut self, point: [f32; 3]) {
        self.0.extend(point);
    }

    fn union(&mut self, other: PyAabb) {
        self.0.union(&other.0);
    }

    fn within(&self, other: PyAabb, distance: f32) -> bool {
        self.0.within(&other.0, distance)
    }
}

#[pyclass(name = "SymbolId", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct PySymbolId(pub(crate) pdbiox::SymbolId);

#[pymethods]
impl PySymbolId {
    #[new]
    fn new(value: u32) -> Self {
        Self(pdbiox::SymbolId::from_raw(value))
    }

    #[staticmethod]
    fn from_raw(value: u32) -> Self {
        Self(pdbiox::SymbolId::from_raw(value))
    }

    #[getter]
    fn value(&self) -> u32 {
        self.0.get()
    }

    fn get(&self) -> u32 {
        self.0.get()
    }

    fn is_canonical(&self) -> bool {
        self.0.is_canonical()
    }

    fn __int__(&self) -> u32 {
        self.0.get()
    }

    fn __index__(&self) -> u32 {
        self.0.get()
    }

    fn __repr__(&self) -> String {
        format!("SymbolId({})", self.0.get())
    }
}

#[pyclass(name = "AltId", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct PyAltId(pub(crate) pdbiox::AltId);

#[pymethods]
impl PyAltId {
    #[new]
    fn new(value: u32) -> Self {
        Self(pdbiox::AltId::from_raw(value))
    }

    #[staticmethod]
    fn from_raw(value: u32) -> Self {
        Self(pdbiox::AltId::from_raw(value))
    }

    #[classattr]
    #[pyo3(name = "BLANK")]
    fn blank() -> Self {
        Self(pdbiox::AltId::BLANK)
    }

    #[staticmethod]
    fn labelled(symbol: PySymbolId) -> Option<Self> {
        pdbiox::AltId::labelled(symbol.0).map(Self)
    }

    #[getter]
    fn value(&self) -> u32 {
        self.0.get()
    }

    fn get(&self) -> u32 {
        self.0.get()
    }

    fn is_blank(&self) -> bool {
        self.0.is_blank()
    }

    fn symbol(&self) -> Option<PySymbolId> {
        self.0.symbol().map(PySymbolId)
    }

    fn __int__(&self) -> u32 {
        self.0.get()
    }
}

#[pyclass(name = "OptionalI32", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyOptionalI32(pub(crate) pdbiox::core::OptionalI32);

#[pymethods]
impl PyOptionalI32 {
    #[new]
    #[pyo3(signature = (value=None))]
    pub(crate) fn new(value: Option<i32>) -> Self {
        Self(value.map_or(
            pdbiox::core::OptionalI32::NONE,
            pdbiox::core::OptionalI32::some,
        ))
    }

    #[staticmethod]
    fn some(value: i32) -> Self {
        Self(pdbiox::core::OptionalI32::some(value))
    }

    #[classattr]
    #[pyo3(name = "NONE")]
    fn none() -> Self {
        Self(pdbiox::core::OptionalI32::NONE)
    }

    fn get(&self) -> Option<i32> {
        self.0.get()
    }

    fn is_none(&self) -> bool {
        self.0.is_none()
    }

    fn is_some(&self) -> bool {
        self.0.is_some()
    }
}

#[pyclass(name = "OptionalSymbol", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyOptionalSymbol(pub(crate) pdbiox::core::OptionalSymbol);

#[pymethods]
impl PyOptionalSymbol {
    #[new]
    #[pyo3(signature = (value=None))]
    pub(crate) fn new(value: Option<PySymbolId>) -> Self {
        Self(value.map_or(pdbiox::core::OptionalSymbol::NONE, |value| {
            pdbiox::core::OptionalSymbol::some(value.0)
        }))
    }

    #[staticmethod]
    fn some(value: PySymbolId) -> Self {
        Self(pdbiox::core::OptionalSymbol::some(value.0))
    }

    #[classattr]
    #[pyo3(name = "NONE")]
    fn none() -> Self {
        Self(pdbiox::core::OptionalSymbol::NONE)
    }

    fn get(&self) -> Option<PySymbolId> {
        self.0.get().map(PySymbolId)
    }

    fn is_none(&self) -> bool {
        self.0.is_none()
    }

    fn is_some(&self) -> bool {
        self.0.is_some()
    }
}

#[pyclass(name = "Interner", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyInterner(pub(crate) pdbiox::core::Interner);

#[pymethods]
impl PyInterner {
    #[new]
    fn new() -> Self {
        Self(pdbiox::core::Interner::new())
    }

    #[classattr]
    #[pyo3(name = "DEFAULT_LIMIT")]
    fn default_limit() -> u32 {
        pdbiox::core::Interner::DEFAULT_LIMIT
    }

    #[staticmethod]
    fn with_capacity(identifiers: usize) -> Self {
        Self(pdbiox::core::Interner::with_capacity(identifiers))
    }

    fn with_limit(&self, limit: u32) -> Self {
        Self(self.0.clone().with_limit(limit))
    }

    fn intern(&mut self, text: &str) -> PyResult<PySymbolId> {
        self.0
            .intern(text)
            .map(PySymbolId)
            .map_err(|error| DictionaryFull::new_err(error.to_string()))
    }

    fn get(&self, text: &str) -> Option<PySymbolId> {
        self.0.get(text).map(PySymbolId)
    }

    fn resolve(&self, symbol: PySymbolId) -> Option<String> {
        self.0.resolve(symbol.0).map(str::to_owned)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn arena_len(&self) -> usize {
        self.0.arena_len()
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> Vec<(PySymbolId, String)> {
        self.0
            .iter()
            .map(|(symbol, text)| (PySymbolId(symbol), text.to_owned()))
            .collect()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPosition>()?;
    module.add_class::<PyByteSpan>()?;
    module.add_class::<PyCoordinateGeneration>()?;
    module.add_class::<PyAabb>()?;
    module.add_class::<PySymbolId>()?;
    module.add_class::<PyAltId>()?;
    module.add_class::<PyOptionalI32>()?;
    module.add_class::<PyOptionalSymbol>()?;
    module.add_class::<PyInterner>()?;
    module.add("DictionaryFull", module.py().get_type::<DictionaryFull>())?;
    Ok(())
}
