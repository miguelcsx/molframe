//! Reflection tables: reading and writing MTZ, and the columns of reciprocal-space data.

use crate::crystal::PyUnitCell;
use molframe::crystal::{
    ReflectionColumn, ReflectionColumnType, ReflectionTable, ReflectionValue, read_mtz, write_mtz,
};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

/// Reflections: Miller indices and the named columns measured at them.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ReflectionTable",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyReflectionTable {
    inner: Arc<ReflectionTable>,
}

fn real(value: &ReflectionValue) -> f64 {
    match value.as_f64() {
        Some(number) => number,
        None => f64::NAN,
    }
}

#[pymethods]
impl PyReflectionTable {
    /// Builds a table from Miller indices `(n, 3)` and columns, each `(type, values)` with
    /// `values` of length `n` (`nan` marks a value that was not recorded).
    ///
    /// The types are `amplitude`, `intensity`, `standard_deviation`, `phase`, `flag` and
    /// `real`. The Miller index columns are written as `H`, `K` and `L`.
    #[new]
    #[pyo3(signature = (
        hkl,
        columns,
        *,
        cell=None,
        space_group_number=None,
        space_group_name=None,
        title="",
    ))]
    #[allow(clippy::needless_pass_by_value)]
    fn new(
        hkl: &Bound<'_, PyArray2<i32>>,
        columns: BTreeMap<String, (String, Vec<f64>)>,
        cell: Option<&PyUnitCell>,
        space_group_number: Option<i32>,
        space_group_name: Option<String>,
        title: &str,
    ) -> PyResult<Self> {
        let readonly = hkl.readonly();
        let view = readonly.as_array();
        if view.ncols() != 3 {
            return Err(crate::error::value("hkl must have three columns"));
        }
        let mut table = ReflectionTable {
            title: title.into(),
            cell: cell.map(PyUnitCell::parameters),
            space_group_number,
            space_group_name: space_group_name.map(String::into_boxed_str),
            ..ReflectionTable::default()
        };
        for (position, label) in ["H", "K", "L"].into_iter().enumerate() {
            table.columns.push(ReflectionColumn {
                label: label.into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: view
                    .column(position)
                    .iter()
                    .map(|&index| ReflectionValue::Integer(i64::from(index)))
                    .collect(),
                dataset_id: 0,
                mtz_type: Some('H'),
            });
        }
        for (label, (kind, values)) in columns {
            if values.len() != view.nrows() {
                return Err(crate::error::value(format!(
                    "column {label:?} has {} values for {} reflections",
                    values.len(),
                    view.nrows()
                )));
            }
            let column_type: ReflectionColumnType = kind.parse().map_err(crate::error::kernel)?;
            table.columns.push(ReflectionColumn {
                label: label.into_boxed_str(),
                column_type,
                values: values
                    .into_iter()
                    .map(|value| {
                        if value.is_nan() {
                            ReflectionValue::Missing
                        } else {
                            ReflectionValue::Real(value)
                        }
                    })
                    .collect(),
                dataset_id: 0,
                mtz_type: None,
            });
        }
        table
            .resolve_symmetry_operations()
            .map_err(crate::error::kernel)?;
        table.validate().map_err(crate::error::kernel)?;
        Ok(Self {
            inner: Arc::new(table),
        })
    }

    /// The title.
    #[getter]
    fn title(&self) -> String {
        self.inner.title.to_string()
    }

    /// The unit cell, when the file gave one.
    #[getter]
    fn cell(&self) -> PyResult<Option<PyUnitCell>> {
        self.inner
            .cell
            .as_ref()
            .map(PyUnitCell::from_cell)
            .transpose()
    }

    /// The International Tables space-group number, when given.
    #[getter]
    fn space_group_number(&self) -> Option<i32> {
        self.inner.space_group_number
    }

    /// The Hermann–Mauguin name, when given.
    #[getter]
    fn space_group_name(&self) -> Option<String> {
        self.inner
            .space_group_name
            .as_ref()
            .map(ToString::to_string)
    }

    /// The number of reflections.
    #[getter]
    fn row_count(&self) -> usize {
        self.inner.row_count()
    }

    fn __len__(&self) -> usize {
        self.inner.row_count()
    }

    /// The column labels, in file order.
    #[getter]
    fn labels(&self) -> Vec<String> {
        self.inner
            .columns
            .iter()
            .map(|column| column.label.to_string())
            .collect()
    }

    /// Each column's type: `amplitude`, `intensity`, `standard_deviation`, `phase`, `flag`,
    /// `real`, `miller_index` or `text`.
    #[getter]
    fn types(&self) -> BTreeMap<String, &'static str> {
        self.inner
            .columns
            .iter()
            .map(|column| (column.label.to_string(), column.column_type.name()))
            .collect()
    }

    /// The Miller indices, `(n, 3)`.
    #[getter]
    fn miller_indices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i32>>> {
        let rows = self.inner.miller_indices().map_err(crate::error::kernel)?;
        let mut flat = Vec::with_capacity(rows.len() * 3);
        let count = rows.len();
        for row in rows {
            flat.extend(row.map_err(crate::error::kernel)?);
        }
        flat.into_pyarray(py).reshape((count, 3))
    }

    /// One column as floats, `nan` where the value was not recorded; the text columns have
    /// no numeric form and are refused.
    fn column<'py>(&self, py: Python<'py>, label: &str) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let column = self
            .inner
            .columns
            .iter()
            .find(|column| column.label.eq_ignore_ascii_case(label))
            .ok_or_else(|| crate::error::key(label))?;
        if column.column_type == ReflectionColumnType::Text {
            return Err(crate::error::value(format!("column {label:?} is text")));
        }
        let values: Vec<f64> = column.values.iter().map(real).collect();
        Ok(values.into_pyarray(py))
    }

    /// The history lines the file carries, in order.
    #[getter]
    fn history(&self) -> Vec<String> {
        self.inner.history.iter().map(ToString::to_string).collect()
    }

    /// The table as MTZ bytes.
    fn to_mtz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let table = Arc::clone(&self.inner);
        let bytes = py
            .detach(move || write_mtz(&table))
            .map_err(crate::error::kernel)?;
        Ok(PyBytes::new(py, &bytes))
    }

    fn __repr__(&self) -> String {
        format!(
            "ReflectionTable(reflections={}, columns={})",
            self.inner.row_count(),
            self.inner.columns.len()
        )
    }
}

/// Reads an MTZ file.
#[pyfunction]
#[pyo3(name = "read_mtz")]
fn read_mtz_file(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<PyReflectionTable> {
    let bytes: Vec<u8> = if let Ok(path) = source.extract::<PathBuf>() {
        std::fs::read(&path).map_err(|error| {
            crate::error::kernel(
                molframe::Diagnostic::new(molframe::Code::E7101)
                    .with_message(format!("{}: {error}", path.display())),
            )
        })?
    } else {
        source
            .extract::<pyo3::pybacked::PyBackedBytes>()
            .map_err(|_| crate::error::type_error("source must be a path or bytes"))?
            .to_vec()
    };
    let table = py
        .detach(move || read_mtz(&bytes))
        .map_err(crate::error::kernel)?;
    Ok(PyReflectionTable {
        inner: Arc::new(table),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyReflectionTable>()?;
    module.add_function(wrap_pyfunction!(read_mtz_file, module)?)
}
