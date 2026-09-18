//! Typed Python projection of the generic native encoded-column storage.

use crate::core_values::PySymbolId;
use numpy::ndarray::ArrayView1;
use numpy::{Element, PyArray1, PyArrayMethods};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList};

#[pyclass(name = "ColumnKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyColumnKind {
    U8,
    U16,
    U32,
    U64,
    I8,
    I16,
    I32,
    I64,
    F32,
    SymbolId,
}

#[derive(Clone)]
enum ColumnData {
    U8(molframe::EncodedColumn<u8>),
    U16(molframe::EncodedColumn<u16>),
    U32(molframe::EncodedColumn<u32>),
    U64(molframe::EncodedColumn<u64>),
    I8(molframe::EncodedColumn<i8>),
    I16(molframe::EncodedColumn<i16>),
    I32(molframe::EncodedColumn<i32>),
    I64(molframe::EncodedColumn<i64>),
    F32(molframe::EncodedColumn<f32>),
    SymbolId(molframe::EncodedColumn<molframe::SymbolId>),
}

#[pyclass(name = "EncodedColumn", from_py_object)]
#[derive(Clone)]
pub(crate) struct PyEncodedColumn {
    kind: PyColumnKind,
    data: ColumnData,
}

#[pymethods]
impl PyEncodedColumn {
    #[staticmethod]
    fn plain(values: &Bound<'_, PyAny>, kind: PyColumnKind) -> PyResult<Self> {
        Self::build(values, kind, false)
    }

    #[staticmethod]
    fn encode(values: &Bound<'_, PyAny>, kind: PyColumnKind) -> PyResult<Self> {
        Self::build(values, kind, true)
    }

    #[getter]
    fn kind(&self) -> PyColumnKind {
        self.kind
    }

    fn len(&self) -> usize {
        self.data.len()
    }

    fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    fn is_constant(&self) -> bool {
        self.data.is_constant()
    }

    fn get(&self, py: Python<'_>, position: u32) -> PyResult<Option<Py<PyAny>>> {
        match &self.data {
            ColumnData::U8(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::U16(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::U32(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::U64(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::I8(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::I16(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::I32(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::I64(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::F32(value) => value.get(position).map(|v| v.into_py_any(py)).transpose(),
            ColumnData::SymbolId(value) => value
                .get(position)
                .map(|v| Py::new(py, PySymbolId(v)).map(Py::into_any))
                .transpose(),
        }
    }

    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let values = PyList::empty(py);
        match &self.data {
            ColumnData::U8(value) => append_values(py, &values, value.iter())?,
            ColumnData::U16(value) => append_values(py, &values, value.iter())?,
            ColumnData::U32(value) => append_values(py, &values, value.iter())?,
            ColumnData::U64(value) => append_values(py, &values, value.iter())?,
            ColumnData::I8(value) => append_values(py, &values, value.iter())?,
            ColumnData::I16(value) => append_values(py, &values, value.iter())?,
            ColumnData::I32(value) => append_values(py, &values, value.iter())?,
            ColumnData::I64(value) => append_values(py, &values, value.iter())?,
            ColumnData::F32(value) => append_values(py, &values, value.iter())?,
            ColumnData::SymbolId(value) => {
                for item in value {
                    values.append(Py::new(py, PySymbolId(item))?)?;
                }
            }
        }
        Ok(values)
    }

    #[pyo3(name = "iter")]
    fn iter_values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        self.to_list(py)
    }

    /// Returns a read-only `NumPy` view only for a physically plain numeric
    /// column. Encoded columns and symbol columns return `None` rather than
    /// materialising a misleading copy behind an `as_slice` name.
    fn as_slice<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        let owner = Bound::new(py, self.clone())?;
        let plain = {
            let value = owner.borrow();
            match &value.data {
                ColumnData::U8(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::U8(values.as_ptr(), values.len())),
                ColumnData::U16(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::U16(values.as_ptr(), values.len())),
                ColumnData::U32(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::U32(values.as_ptr(), values.len())),
                ColumnData::U64(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::U64(values.as_ptr(), values.len())),
                ColumnData::I8(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::I8(values.as_ptr(), values.len())),
                ColumnData::I16(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::I16(values.as_ptr(), values.len())),
                ColumnData::I32(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::I32(values.as_ptr(), values.len())),
                ColumnData::I64(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::I64(values.as_ptr(), values.len())),
                ColumnData::F32(column) => column
                    .as_slice()
                    .map(|values| PlainSlice::F32(values.as_ptr(), values.len())),
                ColumnData::SymbolId(_) => None,
            }
        };
        let Some(plain) = plain else {
            return Ok(None);
        };
        let owner = owner.into_any();
        Ok(Some(match plain {
            PlainSlice::U8(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::U16(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::U32(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::U64(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::I8(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::I16(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::I32(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::I64(pointer, length) => plain_array(owner, pointer, length),
            PlainSlice::F32(pointer, length) => plain_array(owner, pointer, length),
        }))
    }

    fn __len__(&self) -> usize {
        self.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "EncodedColumn(kind={:?}, len={}, constant={})",
            self.kind,
            self.len(),
            self.is_constant()
        )
    }
}

enum PlainSlice {
    U8(*const u8, usize),
    U16(*const u16, usize),
    U32(*const u32, usize),
    U64(*const u64, usize),
    I8(*const i8, usize),
    I16(*const i16, usize),
    I32(*const i32, usize),
    I64(*const i64, usize),
    F32(*const f32, usize),
}

fn plain_array<T: Element>(
    owner: Bound<'_, PyAny>,
    pointer: *const T,
    length: usize,
) -> Bound<'_, PyAny> {
    // SAFETY: `pointer` and `length` were read from the cloned column retained
    // by `owner`; no mutable path exists on this frozen view.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: the NumPy base retains the cloned Rust column for the view life.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner) };
    array.readwrite().make_nonwriteable();
    array.into_any()
}

impl PyEncodedColumn {
    fn build(values: &Bound<'_, PyAny>, kind: PyColumnKind, encoded: bool) -> PyResult<Self> {
        macro_rules! build_column {
            ($type:ty, $variant:ident) => {{
                let values = values.extract::<Vec<$type>>()?;
                let column = if encoded {
                    molframe::EncodedColumn::encode(&values)
                } else {
                    molframe::EncodedColumn::plain(&values)
                };
                Ok(Self {
                    kind,
                    data: ColumnData::$variant(column),
                })
            }};
        }

        match kind {
            PyColumnKind::U8 => build_column!(u8, U8),
            PyColumnKind::U16 => build_column!(u16, U16),
            PyColumnKind::U32 => build_column!(u32, U32),
            PyColumnKind::U64 => build_column!(u64, U64),
            PyColumnKind::I8 => build_column!(i8, I8),
            PyColumnKind::I16 => build_column!(i16, I16),
            PyColumnKind::I32 => build_column!(i32, I32),
            PyColumnKind::I64 => build_column!(i64, I64),
            PyColumnKind::F32 => build_column!(f32, F32),
            PyColumnKind::SymbolId => {
                let values = values.extract::<Vec<PySymbolId>>()?;
                let values = values.into_iter().map(|value| value.0).collect::<Vec<_>>();
                let column = if encoded {
                    molframe::EncodedColumn::encode(&values)
                } else {
                    molframe::EncodedColumn::plain(&values)
                };
                Ok(Self {
                    kind,
                    data: ColumnData::SymbolId(column),
                })
            }
        }
    }
}

trait ColumnLength {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
    fn is_constant(&self) -> bool;
}

impl<T: molframe::core::column::ColumnValue> ColumnLength for molframe::EncodedColumn<T> {
    fn len(&self) -> usize {
        molframe::EncodedColumn::len(self)
    }

    fn is_empty(&self) -> bool {
        molframe::EncodedColumn::is_empty(self)
    }

    fn is_constant(&self) -> bool {
        molframe::EncodedColumn::is_constant(self)
    }
}

impl ColumnLength for ColumnData {
    fn len(&self) -> usize {
        match self {
            Self::U8(value) => value.len(),
            Self::U16(value) => value.len(),
            Self::U32(value) => value.len(),
            Self::U64(value) => value.len(),
            Self::I8(value) => value.len(),
            Self::I16(value) => value.len(),
            Self::I32(value) => value.len(),
            Self::I64(value) => value.len(),
            Self::F32(value) => value.len(),
            Self::SymbolId(value) => value.len(),
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn is_constant(&self) -> bool {
        match self {
            Self::U8(value) => value.is_constant(),
            Self::U16(value) => value.is_constant(),
            Self::U32(value) => value.is_constant(),
            Self::U64(value) => value.is_constant(),
            Self::I8(value) => value.is_constant(),
            Self::I16(value) => value.is_constant(),
            Self::I32(value) => value.is_constant(),
            Self::I64(value) => value.is_constant(),
            Self::F32(value) => value.is_constant(),
            Self::SymbolId(value) => value.is_constant(),
        }
    }
}

fn append_values<T>(
    py: Python<'_>,
    target: &Bound<'_, PyList>,
    values: impl Iterator<Item = T>,
) -> PyResult<()>
where
    T: for<'a> IntoPyObject<'a>,
{
    for value in values {
        target.append(value.into_py_any(py)?)?;
    }
    Ok(())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyColumnKind>()?;
    module.add_class::<PyEncodedColumn>()?;
    Ok(())
}
