//! Typed, structure-local atom annotations.

use crate::core_columns::PyPresence;
use crate::core_values::PySymbolId;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyFloat, PyInt, PyList};

#[derive(Clone, Debug)]
enum ColumnInner {
    Boolean(molframe::core::AnnotationColumn<bool>),
    Integer(molframe::core::AnnotationColumn<i64>),
    Real(molframe::core::AnnotationColumn<f64>),
    Symbol(molframe::core::AnnotationColumn<molframe::core::SymbolId>),
}

#[pyclass(name = "AnnotationColumn", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAnnotationColumn {
    inner: ColumnInner,
}

#[derive(Clone, Debug)]
enum AnnotationInner {
    Boolean(molframe::core::AnnotationColumn<bool>),
    Integer(molframe::core::AnnotationColumn<i64>),
    Real(molframe::core::AnnotationColumn<f64>),
    Symbol(molframe::core::AnnotationColumn<molframe::core::SymbolId>),
}

#[pyclass(name = "AtomAnnotation", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAtomAnnotation {
    inner: AnnotationInner,
}

#[pyclass(name = "AtomAnnotations", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAtomAnnotations(pub(crate) molframe::core::AtomAnnotations);

fn capacity_error(error: molframe::core::CapacityError) -> PyErr {
    pyo3::exceptions::PyOverflowError::new_err(error.to_string())
}

fn entries<T: Copy>(
    values: Vec<T>,
    presences: Vec<PyPresence>,
) -> PyResult<molframe::core::AnnotationColumn<T>> {
    if values.len() != presences.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "values and presences must have equal length",
        ));
    }
    molframe::core::AnnotationColumn::from_entries(
        values
            .into_iter()
            .zip(presences.into_iter().map(Into::into)),
    )
    .map_err(capacity_error)
}

#[pymethods]
impl PyAnnotationColumn {
    #[staticmethod]
    fn boolean(values: Vec<bool>) -> PyResult<Self> {
        molframe::core::AnnotationColumn::from_values(values)
            .map(|inner| Self {
                inner: ColumnInner::Boolean(inner),
            })
            .map_err(capacity_error)
    }

    #[staticmethod]
    fn integer(values: Vec<i64>) -> PyResult<Self> {
        molframe::core::AnnotationColumn::from_values(values)
            .map(|inner| Self {
                inner: ColumnInner::Integer(inner),
            })
            .map_err(capacity_error)
    }

    #[staticmethod]
    fn real(values: Vec<f64>) -> PyResult<Self> {
        molframe::core::AnnotationColumn::from_values(values)
            .map(|inner| Self {
                inner: ColumnInner::Real(inner),
            })
            .map_err(capacity_error)
    }

    #[staticmethod]
    fn symbol(values: Vec<PySymbolId>) -> PyResult<Self> {
        let values = values.into_iter().map(|value| value.0).collect();
        molframe::core::AnnotationColumn::from_values(values)
            .map(|inner| Self {
                inner: ColumnInner::Symbol(inner),
            })
            .map_err(capacity_error)
    }

    #[staticmethod]
    fn boolean_entries(values: Vec<bool>, presences: Vec<PyPresence>) -> PyResult<Self> {
        entries(values, presences).map(|inner| Self {
            inner: ColumnInner::Boolean(inner),
        })
    }

    #[staticmethod]
    fn integer_entries(values: Vec<i64>, presences: Vec<PyPresence>) -> PyResult<Self> {
        entries(values, presences).map(|inner| Self {
            inner: ColumnInner::Integer(inner),
        })
    }

    #[staticmethod]
    fn real_entries(values: Vec<f64>, presences: Vec<PyPresence>) -> PyResult<Self> {
        entries(values, presences).map(|inner| Self {
            inner: ColumnInner::Real(inner),
        })
    }

    #[staticmethod]
    fn symbol_entries(values: Vec<PySymbolId>, presences: Vec<PyPresence>) -> PyResult<Self> {
        let values = values.into_iter().map(|value| value.0).collect();
        entries(values, presences).map(|inner| Self {
            inner: ColumnInner::Symbol(inner),
        })
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.inner {
            ColumnInner::Boolean(_) => "boolean",
            ColumnInner::Integer(_) => "integer",
            ColumnInner::Real(_) => "real",
            ColumnInner::Symbol(_) => "symbol",
        }
    }

    #[getter]
    fn len(&self) -> u32 {
        match &self.inner {
            ColumnInner::Boolean(column) => column.len(),
            ColumnInner::Integer(column) => column.len(),
            ColumnInner::Real(column) => column.len(),
            ColumnInner::Symbol(column) => column.len(),
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn presence(&self, atom: u32) -> PyPresence {
        match &self.inner {
            ColumnInner::Boolean(column) => column.presence(atom),
            ColumnInner::Integer(column) => column.presence(atom),
            ColumnInner::Real(column) => column.presence(atom),
            ColumnInner::Symbol(column) => column.presence(atom),
        }
        .into()
    }

    fn get(&self, py: Python<'_>, atom: u32) -> PyResult<Option<(Py<PyAny>, PyPresence)>> {
        let value = match &self.inner {
            ColumnInner::Boolean(column) => column.get(atom).map(|(value, presence)| {
                (
                    PyBool::new(py, value).to_owned().unbind().into_any(),
                    presence.into(),
                )
            }),
            ColumnInner::Integer(column) => column.get(atom).map(|(value, presence)| {
                (PyInt::new(py, value).unbind().into_any(), presence.into())
            }),
            ColumnInner::Real(column) => column.get(atom).map(|(value, presence)| {
                (PyFloat::new(py, value).unbind().into_any(), presence.into())
            }),
            ColumnInner::Symbol(column) => column
                .get(atom)
                .map(|(value, presence)| {
                    Py::new(py, PySymbolId(value)).map(|value| (value.into_any(), presence.into()))
                })
                .transpose()?,
        };
        Ok(value)
    }

    #[getter]
    fn values(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.inner {
            ColumnInner::Boolean(column) => {
                Ok(PyList::new(py, column.values())?.unbind().into_any())
            }
            ColumnInner::Integer(column) => {
                Ok(PyList::new(py, column.values())?.unbind().into_any())
            }
            ColumnInner::Real(column) => Ok(PyList::new(py, column.values())?.unbind().into_any()),
            ColumnInner::Symbol(column) => {
                let values = column
                    .values()
                    .iter()
                    .copied()
                    .map(|value| Py::new(py, PySymbolId(value)).map(Py::into_any))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok(PyList::new(py, values)?.unbind().into_any())
            }
        }
    }
}

#[pymethods]
impl PyAtomAnnotation {
    #[staticmethod]
    fn boolean(column: PyAnnotationColumn) -> PyResult<Self> {
        match column.inner {
            ColumnInner::Boolean(column) => Ok(Self {
                inner: AnnotationInner::Boolean(column),
            }),
            _ => Err(pyo3::exceptions::PyTypeError::new_err(
                "expected a boolean annotation column",
            )),
        }
    }

    #[staticmethod]
    fn integer(column: PyAnnotationColumn) -> PyResult<Self> {
        match column.inner {
            ColumnInner::Integer(column) => Ok(Self {
                inner: AnnotationInner::Integer(column),
            }),
            _ => Err(pyo3::exceptions::PyTypeError::new_err(
                "expected an integer annotation column",
            )),
        }
    }

    #[staticmethod]
    fn real(column: PyAnnotationColumn) -> PyResult<Self> {
        match column.inner {
            ColumnInner::Real(column) => Ok(Self {
                inner: AnnotationInner::Real(column),
            }),
            _ => Err(pyo3::exceptions::PyTypeError::new_err(
                "expected a real annotation column",
            )),
        }
    }

    #[staticmethod]
    fn symbol(column: PyAnnotationColumn) -> PyResult<Self> {
        match column.inner {
            ColumnInner::Symbol(column) => Ok(Self {
                inner: AnnotationInner::Symbol(column),
            }),
            _ => Err(pyo3::exceptions::PyTypeError::new_err(
                "expected a symbol annotation column",
            )),
        }
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.inner {
            AnnotationInner::Boolean(_) => "boolean",
            AnnotationInner::Integer(_) => "integer",
            AnnotationInner::Real(_) => "real",
            AnnotationInner::Symbol(_) => "symbol",
        }
    }

    fn len(&self) -> u32 {
        match &self.inner {
            AnnotationInner::Boolean(column) => column.len(),
            AnnotationInner::Integer(column) => column.len(),
            AnnotationInner::Real(column) => column.len(),
            AnnotationInner::Symbol(column) => column.len(),
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn column(&self) -> PyAnnotationColumn {
        let inner = match &self.inner {
            AnnotationInner::Boolean(column) => ColumnInner::Boolean(column.clone()),
            AnnotationInner::Integer(column) => ColumnInner::Integer(column.clone()),
            AnnotationInner::Real(column) => ColumnInner::Real(column.clone()),
            AnnotationInner::Symbol(column) => ColumnInner::Symbol(column.clone()),
        };
        PyAnnotationColumn { inner }
    }
}

#[pymethods]
impl PyAtomAnnotations {
    #[new]
    fn new() -> Self {
        Self(molframe::core::AtomAnnotations::default())
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn get(&self, name: &str) -> PyResult<Option<PyAtomAnnotation>> {
        self.0
            .get(name)
            .cloned()
            .map(to_python_annotation)
            .transpose()
    }

    fn insert(
        &mut self,
        name: &str,
        annotation: PyAtomAnnotation,
    ) -> PyResult<Option<PyAtomAnnotation>> {
        self.0
            .insert(name, annotation.into())
            .map(to_python_annotation)
            .transpose()
    }

    fn remove(&mut self, name: &str) -> PyResult<Option<PyAtomAnnotation>> {
        self.0.remove(name).map(to_python_annotation).transpose()
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> PyResult<Vec<(String, PyAtomAnnotation)>> {
        self.0
            .iter()
            .map(|(name, value)| {
                to_python_annotation(value.clone()).map(|value| (name.to_owned(), value))
            })
            .collect()
    }
}

fn to_python_annotation(value: molframe::core::AtomAnnotation) -> PyResult<PyAtomAnnotation> {
    let inner = match value {
        molframe::core::AtomAnnotation::Boolean(column) => AnnotationInner::Boolean(column),
        molframe::core::AtomAnnotation::Integer(column) => AnnotationInner::Integer(column),
        molframe::core::AtomAnnotation::Real(column) => AnnotationInner::Real(column),
        molframe::core::AtomAnnotation::Symbol(column) => AnnotationInner::Symbol(column),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "unsupported non-exhaustive atom annotation variant",
            ));
        }
    };
    Ok(PyAtomAnnotation { inner })
}

impl From<PyAtomAnnotation> for molframe::core::AtomAnnotation {
    fn from(value: PyAtomAnnotation) -> Self {
        match value.inner {
            AnnotationInner::Boolean(column) => Self::Boolean(column),
            AnnotationInner::Integer(column) => Self::Integer(column),
            AnnotationInner::Real(column) => Self::Real(column),
            AnnotationInner::Symbol(column) => Self::Symbol(column),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAnnotationColumn>()?;
    module.add_class::<PyAtomAnnotation>()?;
    module.add_class::<PyAtomAnnotations>()?;
    module.add("SEGMENT_ID_ANNOTATION", molframe::SEGMENT_ID_ANNOTATION)?;
    module.add(
        "PARTIAL_CHARGE_ANNOTATION",
        molframe::PARTIAL_CHARGE_ANNOTATION,
    )?;
    module.add("ATOM_RADIUS_ANNOTATION", molframe::ATOM_RADIUS_ANNOTATION)?;
    module.add(
        "AUTODOCK_TYPE_ANNOTATION",
        molframe::AUTODOCK_TYPE_ANNOTATION,
    )?;
    module.add(
        "COMPONENT_KIND_ANNOTATION",
        molframe::COMPONENT_KIND_ANNOTATION,
    )?;
    module.add(
        "POLYMER_ATOM_ROLE_ANNOTATION",
        molframe::POLYMER_ATOM_ROLE_ANNOTATION,
    )?;
    module.add(
        "AROMATIC_ATOM_ANNOTATION",
        molframe::AROMATIC_ATOM_ANNOTATION,
    )?;
    module.add(
        "FORMAL_CHARGE_ANNOTATION",
        molframe::FORMAL_CHARGE_ANNOTATION,
    )?;
    module.add("HBOND_DONOR_ANNOTATION", molframe::HBOND_DONOR_ANNOTATION)?;
    module.add(
        "HBOND_ACCEPTOR_ANNOTATION",
        molframe::HBOND_ACCEPTOR_ANNOTATION,
    )?;
    module.add(
        "STEREO_CONFIGURATION_ANNOTATION",
        molframe::STEREO_CONFIGURATION_ANNOTATION,
    )?;
    module.add("PLDDT_ANNOTATION", molframe::PLDDT_ANNOTATION)?;
    module.add("PAE_ANNOTATION", molframe::PAE_ANNOTATION)?;
    Ok(())
}
