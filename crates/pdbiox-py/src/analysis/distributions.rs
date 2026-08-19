//! Explicit site distributions and membrane components.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::query::PySelection;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

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

impl PyRadialOptions {
    pub(crate) fn native(&self) -> pdbiox::analysis::RadialDistributionOptions {
        rust_options(*self)
    }

    pub(crate) fn from_native(options: pdbiox::analysis::RadialDistributionOptions) -> Self {
        Self {
            minimum: options.minimum_distance,
            maximum: options.maximum_distance,
            bins: options.bins,
            volume: options.volume,
            backend: options.backend.into(),
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

#[pyclass(name = "CentreGroup", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCentreGroup(pub(crate) pdbiox::analysis::CentreGroup);

#[pymethods]
impl PyCentreGroup {
    #[new]
    fn new(atoms: &PySelection) -> Self {
        Self(pdbiox::analysis::CentreGroup {
            atoms: atoms.inner.clone(),
        })
    }

    #[getter]
    fn atoms(&self) -> PySelection {
        PySelection {
            inner: self.0.atoms.clone(),
        }
    }
}

#[pyclass(name = "LeafletOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyLeafletOptions {
    connection_distance: f32,
    backend: PySpatialBackend,
}

#[pymethods]
impl PyLeafletOptions {
    #[new]
    #[pyo3(signature = (connection_distance, *, backend=PySpatialBackend::Auto))]
    fn new(connection_distance: f32, backend: PySpatialBackend) -> Self {
        Self {
            connection_distance,
            backend,
        }
    }
}

impl PyLeafletOptions {
    pub(crate) fn native(&self) -> pdbiox::analysis::LeafletOptions {
        pdbiox::analysis::LeafletOptions {
            connection_distance: self.connection_distance,
            backend: self.backend.into(),
        }
    }

    pub(crate) fn from_native(options: pdbiox::analysis::LeafletOptions) -> Self {
        Self {
            connection_distance: options.connection_distance,
            backend: options.backend.into(),
        }
    }
}

#[pyclass(name = "Leaflet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLeaflet {
    #[pyo3(get)]
    sites: Vec<u32>,
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

impl From<pdbiox::analysis::Leaflet> for PyLeaflet {
    fn from(value: pdbiox::analysis::Leaflet) -> Self {
        Self { sites: value.sites }
    }
}

pub(crate) fn radial_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::RadialBin>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        let values = values
            .into_iter()
            .map(|value| Py::new(py, PyRadialBin::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, values)?.unbind().into_any())
    })
}

pub(crate) fn coordination_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<u32>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        Ok(PyList::new(py, values)?.unbind().into_any())
    })
}

pub(crate) fn leaflets_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::Leaflet>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        let values = values
            .into_iter()
            .map(|value| Py::new(py, PyLeaflet::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, values)?.unbind().into_any())
    })
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
#[pyo3(signature = (structure, left, right, options, *, periodic=false))]
pub(crate) fn radial_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    left: &PySelection,
    right: &PySelection,
    options: &PyRadialOptions,
    periodic: bool,
) -> PyResult<Vec<PyRadialBin>> {
    let structure = structure.structure().clone();
    let left = left.inner.clone();
    let right = right.inner.clone();
    let options = rust_options(*options);
    py.detach(move || {
        let periodic_box = periodic_box(&structure, periodic)?;
        pdbiox::analysis::radial_distribution(
            structure.positions(),
            &left,
            &right,
            options,
            periodic_box.as_ref(),
        )
        .map_err(|error| error.to_string())
    })
    .map(|values| values.into_iter().map(PyRadialBin::from).collect())
    .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (structure, masses, left, right, options, *, periodic=false))]
pub(crate) fn centre_of_mass_radial_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    masses: Vec<f64>,
    left: Vec<PyCentreGroup>,
    right: Vec<PyCentreGroup>,
    options: &PyRadialOptions,
    periodic: bool,
) -> PyResult<Vec<PyRadialBin>> {
    let structure = structure.structure().clone();
    let left: Vec<_> = left.into_iter().map(|value| value.0).collect();
    let right: Vec<_> = right.into_iter().map(|value| value.0).collect();
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
    .map(|values| values.into_iter().map(PyRadialBin::from).collect())
    .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (structure, left, right, shell, *, backend=PySpatialBackend::Auto, periodic=false))]
pub(crate) fn coordination_numbers(
    py: Python<'_>,
    structure: &PyStructure,
    left: &PySelection,
    right: &PySelection,
    shell: (f32, f32),
    backend: PySpatialBackend,
    periodic: bool,
) -> PyResult<Vec<u32>> {
    let structure = structure.structure().clone();
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

#[pyfunction]
#[pyo3(signature = (structure, sites, options, *, periodic=false))]
pub(crate) fn identify_leaflets(
    py: Python<'_>,
    structure: &PyStructure,
    sites: &PySelection,
    options: PyLeafletOptions,
    periodic: bool,
) -> PyResult<Vec<PyLeaflet>> {
    let structure = structure.structure().clone();
    let sites = sites.inner.clone();
    py.detach(move || {
        let periodic_box = periodic_box(&structure, periodic)?;
        pdbiox::analysis::identify_leaflets(
            structure.positions(),
            &sites,
            pdbiox::analysis::LeafletOptions {
                connection_distance: options.connection_distance,
                backend: options.backend.into(),
            },
            periodic_box.as_ref(),
        )
        .map_err(|error| error.to_string())
    })
    .map(|values| {
        values
            .into_iter()
            .map(|value| PyLeaflet { sites: value.sites })
            .collect()
    })
    .map_err(PyValueError::new_err)
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
