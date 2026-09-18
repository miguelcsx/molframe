//! Owned Python policies for native versioned empirical reference data.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyclass(name = "ReferenceDistribution", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceDistribution(pub(super) molframe::validate::ReferenceDistribution);

#[pymethods]
impl PyReferenceDistribution {
    #[staticmethod]
    fn histogram(name: &str, edges: Vec<f64>, weights: Vec<f64>) -> PyResult<Self> {
        molframe::validate::ReferenceDistribution::histogram(name, edges, weights)
            .map(Self)
            .map_err(value_error)
    }

    #[staticmethod]
    fn grid(name: &str, x_edges: Vec<f64>, y_edges: Vec<f64>, weights: Vec<f64>) -> PyResult<Self> {
        molframe::validate::ReferenceDistribution::grid(name, x_edges, y_edges, weights)
            .map(Self)
            .map_err(value_error)
    }
}

#[pyclass(name = "ReferenceLibrary", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceLibrary(pub(super) Arc<molframe::validate::ReferenceLibrary>);

#[pymethods]
impl PyReferenceLibrary {
    #[new]
    fn new(id: &str, version: &str, distributions: Vec<PyReferenceDistribution>) -> PyResult<Self> {
        molframe::validate::ReferenceLibrary::new(
            id,
            version,
            distributions.into_iter().map(|value| value.0),
        )
        .map(|library| Self(Arc::new(library)))
        .map_err(value_error)
    }

    #[getter]
    fn id(&self) -> &str {
        self.0.id()
    }

    #[getter]
    fn version(&self) -> &str {
        self.0.version()
    }
}

#[pyclass(name = "RamachandranBasin", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRamachandranBasin(pub(super) molframe::validate::RamachandranBasin);

impl PyRamachandranBasin {
    fn build(region: molframe::validate::RamachandranRegion, distribution: &str) -> PyResult<Self> {
        molframe::validate::RamachandranBasin::new(region, distribution)
            .map(Self)
            .map_err(value_error)
    }
}

#[pymethods]
impl PyRamachandranBasin {
    #[staticmethod]
    fn alpha_right(distribution: &str) -> PyResult<Self> {
        Self::build(
            molframe::validate::RamachandranRegion::AlphaHelixRight,
            distribution,
        )
    }

    #[staticmethod]
    fn beta_sheet(distribution: &str) -> PyResult<Self> {
        Self::build(
            molframe::validate::RamachandranRegion::BetaSheet,
            distribution,
        )
    }

    #[staticmethod]
    fn alpha_left(distribution: &str) -> PyResult<Self> {
        Self::build(
            molframe::validate::RamachandranRegion::AlphaHelixLeft,
            distribution,
        )
    }
}

#[pyclass(name = "RamachandranOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRamachandranOptions {
    pub(super) references: Arc<molframe::validate::ReferenceLibrary>,
    pub(super) basins: Vec<molframe::validate::RamachandranBasin>,
    pub(super) minimum_probability: f64,
}

#[pymethods]
impl PyRamachandranOptions {
    #[new]
    fn new(
        references: &PyReferenceLibrary,
        basins: Vec<PyRamachandranBasin>,
        minimum_probability: f64,
    ) -> PyResult<Self> {
        let basins = basins.into_iter().map(|value| value.0).collect::<Vec<_>>();
        molframe::validate::RamachandranOptions::new(
            &references.0,
            basins.clone(),
            minimum_probability,
        )
        .map_err(value_error)?;
        Ok(Self {
            references: references.0.clone(),
            basins,
            minimum_probability,
        })
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
