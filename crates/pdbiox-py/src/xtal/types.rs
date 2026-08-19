//! Affine, rational, expression and metadata projections.

use crate::crystallography::PySymmetryOperation;
use crate::geometry::PyRigid;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "AffineTransform", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAffineTransform(pub(crate) pdbiox::xtal::AffineTransform);

#[pymethods]
impl PyAffineTransform {
    #[new]
    fn new(matrix: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self(pdbiox::xtal::AffineTransform::new(matrix, translation))
    }

    #[staticmethod]
    fn identity() -> Self {
        Self(pdbiox::xtal::AffineTransform::IDENTITY)
    }

    #[getter]
    fn matrix(&self) -> [[f64; 3]; 3] {
        self.0.matrix
    }

    #[getter]
    fn translation(&self) -> [f64; 3] {
        self.0.translation
    }

    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    fn apply(&self, position: [f32; 3]) -> [f32; 3] {
        self.0.apply(position)
    }

    fn then(&self, next: &Self) -> Self {
        Self(self.0.then(&next.0))
    }
}

#[pyclass(name = "Rational", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRational {
    #[pyo3(get)]
    numerator: i32,
    #[pyo3(get)]
    denominator: u32,
}

#[pymethods]
impl PyRational {
    #[new]
    fn new(numerator: i32, denominator: u32) -> PyResult<Self> {
        pdbiox::xtal::Rational::new(numerator, denominator)
            .map(Self::from)
            .map_err(value_error)
    }

    #[getter]
    fn value(&self) -> f64 {
        f64::from(self.numerator) / f64::from(self.denominator)
    }
}

#[pyclass(name = "OperExpression", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyOperExpression(pub(crate) pdbiox::xtal::OperExpression);

#[pymethods]
impl PyOperExpression {
    #[staticmethod]
    fn parse(expression: &str) -> PyResult<Self> {
        pdbiox::xtal::OperExpression::parse(expression)
            .map(Self)
            .map_err(value_error)
    }

    #[staticmethod]
    fn parse_with_limit(expression: &str, limit: usize) -> PyResult<Self> {
        pdbiox::xtal::OperExpression::parse_with_limit(expression, limit)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn factors(&self) -> Vec<Vec<String>> {
        self.0
            .factors()
            .iter()
            .map(|factor| factor.iter().map(ToString::to_string).collect())
            .collect()
    }

    #[getter]
    fn combination_count(&self) -> usize {
        self.0.combination_count()
    }
}

#[pyclass(name = "SymmetrySet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySymmetrySet(pub(crate) pdbiox::xtal::SymmetrySet);

#[pymethods]
impl PySymmetrySet {
    #[getter]
    fn hall_number(&self) -> Option<u16> {
        self.0.hall_number
    }
    #[getter]
    fn international_number(&self) -> Option<u16> {
        self.0.international_number
    }
    #[getter]
    fn hermann_mauguin(&self) -> Option<String> {
        self.0.hermann_mauguin.as_deref().map(str::to_owned)
    }
    #[getter]
    fn hall(&self) -> Option<String> {
        self.0.hall.as_deref().map(str::to_owned)
    }
    #[getter]
    fn crystal_system(&self) -> Option<String> {
        self.0.crystal_system.as_deref().map(str::to_owned)
    }
    #[getter]
    fn choice(&self) -> Option<String> {
        self.0.choice.as_deref().map(str::to_owned)
    }
    #[getter]
    fn operations(&self) -> Vec<PySymmetryOperation> {
        self.0
            .operations()
            .iter()
            .map(PySymmetryOperation::from)
            .collect()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[pyclass(name = "Operator", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyOperator {
    #[pyo3(get)]
    pub(crate) id: String,
    #[pyo3(get)]
    pub(crate) transform: PyRigid,
}

#[pyclass(name = "Generator", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGenerator {
    #[pyo3(get)]
    pub(crate) oper_expression: PyOperExpression,
    #[pyo3(get)]
    pub(crate) asym_ids: Vec<String>,
}

#[pyclass(name = "AssemblyDef", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAssemblyDef {
    #[pyo3(get)]
    pub(crate) id: String,
    #[pyo3(get)]
    pub(crate) details: Option<String>,
    #[pyo3(get)]
    pub(crate) method: Option<String>,
    #[pyo3(get)]
    pub(crate) oligomeric: Option<u32>,
    #[pyo3(get)]
    pub(crate) generators: Vec<PyGenerator>,
}

impl From<pdbiox::xtal::Rational> for PyRational {
    fn from(value: pdbiox::xtal::Rational) -> Self {
        Self {
            numerator: value.numerator(),
            denominator: value.denominator(),
        }
    }
}

impl From<pdbiox::xtal::SymmetryOperation> for PySymmetryOperation {
    fn from(value: pdbiox::xtal::SymmetryOperation) -> Self {
        Self::from(&value)
    }
}

impl From<pdbiox::xtal::Operator> for PyOperator {
    fn from(value: pdbiox::xtal::Operator) -> Self {
        Self {
            id: value.id.to_string(),
            transform: PyRigid(value.transform),
        }
    }
}

impl From<pdbiox::xtal::Generator> for PyGenerator {
    fn from(value: pdbiox::xtal::Generator) -> Self {
        Self {
            oper_expression: PyOperExpression(value.oper_expression),
            asym_ids: value.asym_ids.iter().map(ToString::to_string).collect(),
        }
    }
}

impl From<pdbiox::xtal::AssemblyDef> for PyAssemblyDef {
    fn from(value: pdbiox::xtal::AssemblyDef) -> Self {
        Self {
            id: value.id.to_string(),
            details: value.details.map(|value| value.to_string()),
            method: value.method.map(|value| value.to_string()),
            oligomeric: value.oligomeric,
            generators: value.generators.into_iter().map(Into::into).collect(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
