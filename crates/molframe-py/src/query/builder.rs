//! Typed query-builder projections over the shared Rust selection IR.

use super::PyQuery;
use pyo3::prelude::*;

#[pyclass(name = "QueryBuilder", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyQueryBuilder(pub(crate) molframe::query::Builder);

#[pymethods]
impl PyQueryBuilder {
    fn query(&self) -> PyQuery {
        PyQuery {
            inner: molframe::Query::from_builder(self.0.clone()),
        }
    }

    fn __and__(&self, right: &Self) -> Self {
        Self(self.0.clone() & right.0.clone())
    }

    fn __or__(&self, right: &Self) -> Self {
        Self(self.0.clone() | right.0.clone())
    }

    fn __invert__(&self) -> Self {
        Self(!self.0.clone())
    }
}

#[pyclass(name = "ColumnBuilder", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyColumnBuilder(molframe::query::ColumnBuilder);

#[pymethods]
impl PyColumnBuilder {
    fn eq(&self, value: String) -> PyQueryBuilder {
        PyQueryBuilder(self.0.eq(value))
    }

    fn lt(&self, value: f64) -> PyQueryBuilder {
        PyQueryBuilder(self.0.lt(value))
    }

    fn le(&self, value: f64) -> PyQueryBuilder {
        PyQueryBuilder(self.0.le(value))
    }

    fn gt(&self, value: f64) -> PyQueryBuilder {
        PyQueryBuilder(self.0.gt(value))
    }

    fn ge(&self, value: f64) -> PyQueryBuilder {
        PyQueryBuilder(self.0.ge(value))
    }
}

#[pyclass(name = "Columns", frozen)]
pub(crate) struct PyColumns;

#[pymethods]
impl PyColumns {
    #[staticmethod]
    fn all() -> PyQueryBuilder {
        PyQueryBuilder(molframe::query::col::all())
    }

    #[staticmethod]
    fn none() -> PyQueryBuilder {
        PyQueryBuilder(molframe::query::col::none())
    }

    #[staticmethod]
    fn is_protein() -> PyQueryBuilder {
        PyQueryBuilder(molframe::query::col::is_protein())
    }

    #[staticmethod]
    fn name() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::name())
    }

    #[staticmethod]
    fn resname() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::resname())
    }

    #[staticmethod]
    fn chain() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::chain())
    }

    #[staticmethod]
    fn bfactor() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::bfactor())
    }

    #[staticmethod]
    fn occupancy() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::occupancy())
    }

    #[staticmethod]
    fn x() -> PyColumnBuilder {
        PyColumnBuilder(molframe::query::col::x())
    }

    #[staticmethod]
    fn within(radius: f32, target: &PyQueryBuilder) -> PyQueryBuilder {
        PyQueryBuilder(molframe::query::col::within(radius, target.0.clone()))
    }

    #[staticmethod]
    fn by_residue(target: &PyQueryBuilder) -> PyQueryBuilder {
        PyQueryBuilder(molframe::query::col::by_residue(target.0.clone()))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyQueryBuilder>()?;
    module.add_class::<PyColumnBuilder>()?;
    module.add_class::<PyColumns>()?;
    module.add("Builder", module.getattr("QueryBuilder")?)?;
    module.add("col", Py::new(module.py(), PyColumns)?)?;
    Ok(())
}
