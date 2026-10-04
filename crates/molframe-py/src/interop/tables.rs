//! Arrow export for bonds and for named result columns.

use super::stream::arrow_capsule;
use crate::bindings::PyStructure;
use molframe::interop::{BondTable, Column, ColumnTable};
use numpy::PyReadonlyArray1;
use pyo3::prelude::*;
use pyo3::types::PyCapsule;

/// The chemical bonds of a structure, as an Arrow table.
///
/// Columns are `bond_index`, `atom_a`, `atom_b`, `order` and `provenance`; the
/// table is read by any Arrow consumer through `__arrow_c_stream__`.
#[derive(Clone, Debug)]
#[pyclass(name = "BondTable", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyBondTable {
    parent: PyStructure,
}

impl PyBondTable {
    pub(crate) const fn new(parent: PyStructure) -> Self {
        Self { parent }
    }
}

#[pymethods]
impl PyBondTable {
    fn __len__(&self) -> usize {
        self.parent.inner.engine().data().bonds.len()
    }

    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        arrow_capsule(
            py,
            BondTable::new(self.parent.inner.engine()).arrow_stream(),
        )
    }

    fn __repr__(&self) -> String {
        format!("BondTable(bonds={})", self.__len__())
    }
}

/// One column of an analysis result, as the native Arrow column it exports to.
fn column_of(array: &Bound<'_, PyAny>) -> PyResult<Column> {
    macro_rules! carried {
        ($element:ty, $variant:ident) => {
            if let Ok(values) = array.extract::<PyReadonlyArray1<'_, $element>>() {
                return Ok(Column::$variant(values.as_slice()?.to_vec()));
            }
        };
    }
    carried!(u8, U8);
    carried!(u32, U32);
    carried!(i32, I32);
    carried!(i64, I64);
    carried!(f32, F32);
    carried!(f64, F64);
    Err(crate::error::type_error(
        "a column has an element type that Arrow export does not carry",
    ))
}

/// The Arrow stream of named result columns.
pub(crate) fn columns_stream<'py>(
    py: Python<'py>,
    columns: &[(&'static str, Py<PyAny>)],
) -> PyResult<Bound<'py, PyCapsule>> {
    let native = columns
        .iter()
        .map(|(name, array)| Ok(((*name).to_owned(), column_of(array.bind(py))?)))
        .collect::<PyResult<Vec<_>>>()?;
    let table = ColumnTable::new(native).map_err(crate::error::kernel)?;
    arrow_capsule(py, table.arrow_stream())
}
