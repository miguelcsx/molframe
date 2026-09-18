//! Native reciprocal-space table bindings.
//!
//! Reflection values remain owned by Rust.  Heterogeneous columns can be
//! inspected as typed values; dense numeric materialisation is opt-in through
//! `numeric_values`, because an MTZ table cannot always be represented as one
//! zero-copy `NumPy` dtype.

use super::reflection_metadata::PyReflectionMetadata;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyFloat, PyInt, PyString};

#[pyclass(name = "ReflectionColumnType", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyReflectionColumnType {
    MillerIndex,
    Amplitude,
    Intensity,
    StandardDeviation,
    Phase,
    Flag,
    Real,
    Text,
}

#[pyclass(name = "ReflectionValue", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReflectionValue(pub(crate) molframe::xtal::ReflectionValue);

#[pymethods]
impl PyReflectionValue {
    #[new]
    #[pyo3(signature = (kind, value=None))]
    fn new(py: Python<'_>, kind: &str, value: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let value = match kind {
            "missing" => molframe::xtal::ReflectionValue::Missing,
            "inapplicable" => molframe::xtal::ReflectionValue::Inapplicable,
            "integer" => {
                molframe::xtal::ReflectionValue::Integer(required(value, kind)?.extract()?)
            }
            "real" => molframe::xtal::ReflectionValue::Real(required(value, kind)?.extract()?),
            "text" => molframe::xtal::ReflectionValue::Text(
                required(value, kind)?.extract::<String>()?.into(),
            ),
            _ => {
                return Err(PyValueError::new_err(
                    "kind must be missing, inapplicable, integer, real, or text",
                ));
            }
        };
        let _ = py;
        Ok(Self(value))
    }

    #[staticmethod]
    fn missing() -> Self {
        Self(molframe::xtal::ReflectionValue::Missing)
    }

    #[staticmethod]
    fn inapplicable() -> Self {
        Self(molframe::xtal::ReflectionValue::Inapplicable)
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            molframe::xtal::ReflectionValue::Missing => "missing",
            molframe::xtal::ReflectionValue::Inapplicable => "inapplicable",
            molframe::xtal::ReflectionValue::Integer(_) => "integer",
            molframe::xtal::ReflectionValue::Real(_) => "real",
            molframe::xtal::ReflectionValue::Text(_) => "text",
        }
    }

    #[getter]
    fn value(&self, py: Python<'_>) -> Py<PyAny> {
        match &self.0 {
            molframe::xtal::ReflectionValue::Missing
            | molframe::xtal::ReflectionValue::Inapplicable => py.None(),
            molframe::xtal::ReflectionValue::Integer(value) => {
                PyInt::new(py, *value).unbind().into_any()
            }
            molframe::xtal::ReflectionValue::Real(value) => {
                PyFloat::new(py, *value).unbind().into_any()
            }
            molframe::xtal::ReflectionValue::Text(value) => {
                PyString::new(py, value).unbind().into_any()
            }
        }
    }

    fn as_f64(&self) -> Option<f64> {
        self.0.as_f64()
    }

    fn as_i32(&self) -> Option<i32> {
        self.0.as_i32()
    }
}

#[pyclass(name = "ReflectionColumn", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReflectionColumn(pub(crate) molframe::xtal::ReflectionColumn);

#[pymethods]
impl PyReflectionColumn {
    #[new]
    #[pyo3(signature = (label, column_type, values, *, dataset_id=0, mtz_type=None))]
    fn new(
        label: String,
        column_type: PyReflectionColumnType,
        values: &Bound<'_, PyAny>,
        dataset_id: i32,
        mtz_type: Option<String>,
    ) -> PyResult<Self> {
        let mtz_type = mtz_type
            .map(|value| one_character(&value, "mtz_type"))
            .transpose()?;
        Ok(Self(molframe::xtal::ReflectionColumn {
            label: label.into(),
            column_type: column_type.into(),
            values: reflection_values(values)?,
            dataset_id,
            mtz_type,
        }))
    }

    #[getter]
    fn label(&self) -> String {
        self.0.label.to_string()
    }

    #[getter]
    fn column_type(&self) -> PyReflectionColumnType {
        self.0.column_type.into()
    }

    #[getter]
    fn dataset_id(&self) -> i32 {
        self.0.dataset_id
    }

    #[getter]
    fn mtz_type(&self) -> Option<String> {
        self.0.mtz_type.map(|value| value.to_string())
    }

    #[getter]
    fn values(&self) -> Vec<PyReflectionValue> {
        self.0
            .values
            .iter()
            .cloned()
            .map(PyReflectionValue)
            .collect()
    }

    #[pyo3(signature = (*, copy=true))]
    fn numeric_values<'py>(
        &self,
        py: Python<'py>,
        copy: bool,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if !copy {
            return Err(PyValueError::new_err(
                "ReflectionColumn owns a tagged value vector; numeric_values requires copy=True",
            ));
        }
        let values = self
            .0
            .values
            .iter()
            .map(|value| {
                value.as_f64().ok_or_else(|| {
                    PyValueError::new_err("column contains a non-numeric reflection value")
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(values.into_pyarray(py))
    }
}

#[pyclass(name = "ReflectionDataset", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReflectionDataset(pub(crate) molframe::xtal::ReflectionDataset);

#[pymethods]
impl PyReflectionDataset {
    #[new]
    #[pyo3(signature = (id, project, crystal, name, *, wavelength=None, cell=None))]
    fn new(
        id: i32,
        project: String,
        crystal: String,
        name: String,
        wavelength: Option<f64>,
        cell: Option<([f64; 3], [f64; 3])>,
    ) -> Self {
        Self(molframe::xtal::ReflectionDataset {
            id,
            project: project.into(),
            crystal: crystal.into(),
            name: name.into(),
            wavelength,
            cell: cell.map(unit_cell),
        })
    }

    #[getter]
    fn id(&self) -> i32 {
        self.0.id
    }

    #[getter]
    fn project(&self) -> String {
        self.0.project.to_string()
    }

    #[getter]
    fn crystal(&self) -> String {
        self.0.crystal.to_string()
    }

    #[getter]
    fn name(&self) -> String {
        self.0.name.to_string()
    }

    #[getter]
    fn wavelength(&self) -> Option<f64> {
        self.0.wavelength
    }

    #[getter]
    fn cell(&self) -> Option<([f64; 3], [f64; 3])> {
        self.0.cell.map(cell_arrays)
    }
}

#[pyclass(name = "ReflectionTable", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReflectionTable(pub(crate) molframe::xtal::ReflectionTable);

#[pymethods]
impl PyReflectionTable {
    #[new]
    #[pyo3(signature = (columns, *, metadata=None))]
    fn new(
        columns: &Bound<'_, PyAny>,
        metadata: Option<PyRef<'_, PyReflectionMetadata>>,
    ) -> PyResult<Self> {
        let metadata = match metadata {
            Some(metadata) => metadata.clone(),
            None => PyReflectionMetadata::default(),
        };
        Ok(Self(molframe::xtal::ReflectionTable {
            title: metadata.title,
            cell: metadata.cell,
            space_group_number: metadata.space_group_number,
            space_group_name: metadata.space_group_name,
            columns: reflection_columns(columns)?,
            datasets: metadata.datasets,
            history: metadata.history,
            symmetry_operations: metadata.symmetry_operations,
            sort_order: metadata.sort_order,
            resolution_range: metadata.resolution_range,
            missing_value: metadata.missing_value,
            extra_header_records: metadata.extra_header_records,
        }))
    }

    #[getter]
    fn title(&self) -> String {
        self.0.title.to_string()
    }
    #[getter]
    fn cell(&self) -> Option<([f64; 3], [f64; 3])> {
        self.0.cell.map(cell_arrays)
    }
    #[getter]
    fn space_group_number(&self) -> Option<i32> {
        self.0.space_group_number
    }
    #[getter]
    fn space_group_name(&self) -> Option<String> {
        self.0.space_group_name.as_deref().map(str::to_owned)
    }
    #[getter]
    fn columns(&self) -> Vec<PyReflectionColumn> {
        self.0
            .columns
            .iter()
            .cloned()
            .map(PyReflectionColumn)
            .collect()
    }
    #[getter]
    fn datasets(&self) -> Vec<PyReflectionDataset> {
        self.0
            .datasets
            .iter()
            .cloned()
            .map(PyReflectionDataset)
            .collect()
    }
    #[getter]
    fn history(&self) -> Vec<String> {
        self.0.history.iter().map(ToString::to_string).collect()
    }
    #[getter]
    fn symmetry_operations(&self) -> Vec<String> {
        self.0
            .symmetry_operations
            .iter()
            .map(ToString::to_string)
            .collect()
    }
    #[getter]
    fn sort_order(&self) -> [i32; 5] {
        self.0.sort_order
    }
    #[getter]
    fn resolution_range(&self) -> Option<[f64; 2]> {
        self.0.resolution_range
    }
    #[getter]
    fn missing_value(&self) -> Option<f32> {
        self.0.missing_value
    }
    #[getter]
    fn extra_header_records(&self) -> Vec<String> {
        self.0
            .extra_header_records
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    fn row_count(&self) -> usize {
        self.0.row_count()
    }

    fn validate(&self) -> PyResult<()> {
        self.0.validate().map_err(value_error)
    }

    fn column_any(&self, labels: &Bound<'_, PyAny>) -> PyResult<Option<PyReflectionColumn>> {
        let labels: Vec<String> = labels.extract()?;
        let labels = labels.iter().map(String::as_str).collect::<Vec<_>>();
        Ok(self.0.column_any(&labels).cloned().map(PyReflectionColumn))
    }

    fn miller_indices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i32>>> {
        let values = self.0.miller_indices().map_err(value_error)?;
        let rows = values.len();
        let flat = values.into_iter().flatten().collect();
        Array2::from_shape_vec((rows, 3), flat)
            .map(|values| values.into_pyarray(py))
            .map_err(value_error)
    }
}

#[pyfunction]
pub(crate) fn read_mtz(py: Python<'_>, data: &[u8]) -> PyResult<PyReflectionTable> {
    let data = data.to_vec();
    py.detach(move || molframe::xtal::read_mtz(&data))
        .map(PyReflectionTable)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn write_mtz(py: Python<'_>, table: &PyReflectionTable) -> PyResult<Py<PyBytes>> {
    let table = table.0.clone();
    py.detach(move || molframe::xtal::write_mtz(&table))
        .map(|bytes| PyBytes::new(py, &bytes).unbind())
        .map_err(value_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyReflectionMetadata>()?;
    module.add_class::<PyReflectionColumnType>()?;
    module.add_class::<PyReflectionValue>()?;
    module.add_class::<PyReflectionColumn>()?;
    module.add_class::<PyReflectionDataset>()?;
    module.add_class::<PyReflectionTable>()?;
    module.add_function(wrap_pyfunction!(read_mtz, module)?)?;
    module.add_function(wrap_pyfunction!(write_mtz, module)?)?;
    Ok(())
}

impl From<PyReflectionColumnType> for molframe::xtal::ReflectionColumnType {
    fn from(value: PyReflectionColumnType) -> Self {
        match value {
            PyReflectionColumnType::MillerIndex => Self::MillerIndex,
            PyReflectionColumnType::Amplitude => Self::Amplitude,
            PyReflectionColumnType::Intensity => Self::Intensity,
            PyReflectionColumnType::StandardDeviation => Self::StandardDeviation,
            PyReflectionColumnType::Phase => Self::Phase,
            PyReflectionColumnType::Flag => Self::Flag,
            PyReflectionColumnType::Real => Self::Real,
            PyReflectionColumnType::Text => Self::Text,
        }
    }
}

impl From<molframe::xtal::ReflectionColumnType> for PyReflectionColumnType {
    fn from(value: molframe::xtal::ReflectionColumnType) -> Self {
        match value {
            molframe::xtal::ReflectionColumnType::MillerIndex => Self::MillerIndex,
            molframe::xtal::ReflectionColumnType::Amplitude => Self::Amplitude,
            molframe::xtal::ReflectionColumnType::Intensity => Self::Intensity,
            molframe::xtal::ReflectionColumnType::StandardDeviation => Self::StandardDeviation,
            molframe::xtal::ReflectionColumnType::Phase => Self::Phase,
            molframe::xtal::ReflectionColumnType::Flag => Self::Flag,
            molframe::xtal::ReflectionColumnType::Real => Self::Real,
            molframe::xtal::ReflectionColumnType::Text => Self::Text,
        }
    }
}

fn required<'py, 'value>(
    value: Option<&'value Bound<'py, PyAny>>,
    kind: &str,
) -> PyResult<&'value Bound<'py, PyAny>> {
    value.ok_or_else(|| PyValueError::new_err(format!("{kind} reflection value requires value")))
}

fn reflection_values(values: &Bound<'_, PyAny>) -> PyResult<Vec<molframe::xtal::ReflectionValue>> {
    let mut output = Vec::new();
    for value in values.try_iter()? {
        let value = value?;
        let value = value
            .extract::<PyRef<'_, PyReflectionValue>>()
            .map_err(class_error)?;
        output.push(value.0.clone());
    }
    Ok(output)
}

fn reflection_columns(
    values: &Bound<'_, PyAny>,
) -> PyResult<Vec<molframe::xtal::ReflectionColumn>> {
    let mut output = Vec::new();
    for value in values.try_iter()? {
        let value = value?;
        let value = value
            .extract::<PyRef<'_, PyReflectionColumn>>()
            .map_err(class_error)?;
        output.push(value.0.clone());
    }
    Ok(output)
}

fn one_character(value: &str, name: &str) -> PyResult<char> {
    let mut characters = value.chars();
    let character = characters
        .next()
        .ok_or_else(|| PyValueError::new_err(format!("{name} must be one character")))?;
    if characters.next().is_some() {
        return Err(PyValueError::new_err(format!(
            "{name} must be one character"
        )));
    }
    Ok(character)
}

fn unit_cell((lengths, angles): ([f64; 3], [f64; 3])) -> molframe::UnitCell {
    molframe::UnitCell { lengths, angles }
}

fn cell_arrays(cell: molframe::UnitCell) -> ([f64; 3], [f64; 3]) {
    (cell.lengths, cell.angles)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn class_error(error: impl std::fmt::Display) -> PyErr {
    PyTypeError::new_err(error.to_string())
}
