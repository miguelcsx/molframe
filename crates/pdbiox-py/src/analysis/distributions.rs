//! Explicit site distributions and membrane components.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::query::PySelection;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::periodic::periodic_box;

#[pyclass(name = "RadialOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRadialOptions {
    pub(crate) minimum: f32,
    pub(crate) maximum: f32,
    pub(crate) bins: usize,
    pub(crate) volume: f64,
    pub(crate) backend: PySpatialBackend,
}

#[pymethods]
impl PyRadialOptions {
    #[new]
    #[pyo3(signature = (minimum, maximum, bins, volume, *, backend=PySpatialBackend::Auto))]
    fn new(
        minimum: f32,
        maximum: f32,
        bins: usize,
        volume: f64,
        backend: PySpatialBackend,
    ) -> Self {
        Self {
            minimum,
            maximum,
            bins,
            volume,
            backend,
        }
    }
}

#[pyclass(name = "RadialBin", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRadialBin {
    #[pyo3(get)]
    lower: f32,
    #[pyo3(get)]
    upper: f32,
    #[pyo3(get)]
    count: u64,
    #[pyo3(get)]
    distribution: f64,
}

impl From<pdbiox::analysis::RadialBin> for PyRadialBin {
    fn from(bin: pdbiox::analysis::RadialBin) -> Self {
        Self {
            lower: bin.lower,
            upper: bin.upper,
            count: bin.count,
            distribution: bin.distribution,
        }
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (left, right, options, *, periodic=false))]
    fn radial_distribution(
        &self,
        py: Python<'_>,
        left: &PySelection,
        right: &PySelection,
        options: &PyRadialOptions,
        periodic: bool,
    ) -> PyResult<Vec<PyRadialBin>> {
        let structure = self.structure().clone();
        let left = left.inner.clone();
        let right = right.inner.clone();
        let options = *options;
        py.detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            pdbiox::analysis::radial_distribution(
                structure.positions(),
                &left,
                &right,
                pdbiox::analysis::RadialDistributionOptions {
                    minimum_distance: options.minimum,
                    maximum_distance: options.maximum,
                    bins: options.bins,
                    volume: options.volume,
                    backend: options.backend.into(),
                },
                periodic_box.as_ref(),
            )
            .map_err(|error| error.to_string())
        })
        .map(|bins| bins.into_iter().map(Into::into).collect())
        .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (masses, left, right, options, *, periodic=false))]
    fn centre_of_mass_radial_distribution(
        &self,
        py: Python<'_>,
        masses: Vec<f64>,
        left: Vec<PySelection>,
        right: Vec<PySelection>,
        options: &PyRadialOptions,
        periodic: bool,
    ) -> PyResult<Vec<PyRadialBin>> {
        let structure = self.structure().clone();
        let left = groups(left);
        let right = groups(right);
        let options = rust_options(*options);
        py.detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            pdbiox::analysis::centre_of_mass_radial_distribution(
                structure.positions(),
                &masses,
                &left,
                &right,
                options,
                periodic_box.as_ref(),
            )
            .map_err(|error| error.to_string())
        })
        .map(|bins| bins.into_iter().map(Into::into).collect())
        .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (left, right, shell, *, backend=PySpatialBackend::Auto, periodic=false))]
    fn coordination_numbers(
        &self,
        py: Python<'_>,
        left: &PySelection,
        right: &PySelection,
        shell: (f32, f32),
        backend: PySpatialBackend,
        periodic: bool,
    ) -> PyResult<Vec<u32>> {
        let structure = self.structure().clone();
        let left = left.inner.clone();
        let right = right.inner.clone();
        py.detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            pdbiox::analysis::coordination_numbers(
                structure.positions(),
                &left,
                &right,
                shell.0,
                shell.1,
                backend.into(),
                periodic_box.as_ref(),
            )
            .map_err(|error| error.to_string())
        })
        .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (sites, connection_distance, *, backend=PySpatialBackend::Auto, periodic=false))]
    fn leaflets(
        &self,
        py: Python<'_>,
        sites: &PySelection,
        connection_distance: f32,
        backend: PySpatialBackend,
        periodic: bool,
    ) -> PyResult<Vec<Vec<u32>>> {
        let structure = self.structure().clone();
        let sites = sites.inner.clone();
        py.detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            pdbiox::analysis::identify_leaflets(
                structure.positions(),
                &sites,
                pdbiox::analysis::LeafletOptions {
                    connection_distance,
                    backend: backend.into(),
                },
                periodic_box.as_ref(),
            )
            .map(|values| values.into_iter().map(|leaflet| leaflet.sites).collect())
            .map_err(|error| error.to_string())
        })
        .map_err(PyValueError::new_err)
    }
}

#[pyfunction]
pub(crate) fn analyse_centre_of_mass_radial_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    masses: Vec<f64>,
    left: Vec<PySelection>,
    right: Vec<PySelection>,
    options: PyRadialOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let left = groups(left);
    let right = groups(right);
    let options = rust_options(options);
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel = pdbiox::analysis::centre_of_mass_radial_distribution_kernel(
                &masses, &left, &right, options,
            );
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    analysis_with_value(py, analysis, |py, values| {
        let values = values
            .into_iter()
            .map(|value| Py::new(py, PyRadialBin::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(pyo3::types::PyList::new(py, values)?.unbind().into_any())
    })
}

fn groups(values: Vec<PySelection>) -> Vec<pdbiox::analysis::CentreGroup> {
    values
        .into_iter()
        .map(|value| pdbiox::analysis::CentreGroup { atoms: value.inner })
        .collect()
}

fn rust_options(options: PyRadialOptions) -> pdbiox::analysis::RadialDistributionOptions {
    pdbiox::analysis::RadialDistributionOptions {
        minimum_distance: options.minimum,
        maximum_distance: options.maximum,
        bins: options.bins,
        volume: options.volume,
        backend: options.backend.into(),
    }
}
