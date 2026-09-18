//! Explicit site distributions and membrane components.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::core::execution::{PyExecutionContext, default_context};
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
    pub(crate) fn native(&self) -> molframe::analysis::RadialDistributionOptions {
        rust_options(*self)
    }

    pub(crate) fn from_native(options: molframe::analysis::RadialDistributionOptions) -> Self {
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
    pub(super) lower: f32,
    #[pyo3(get)]
    pub(super) upper: f32,
    #[pyo3(get)]
    pub(super) count: u64,
    #[pyo3(get)]
    pub(super) distribution: f64,
}

#[pyclass(name = "CentreGroup", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCentreGroup(pub(crate) molframe::analysis::CentreGroup);

#[pymethods]
impl PyCentreGroup {
    #[new]
    fn new(atoms: &PySelection) -> Self {
        Self(molframe::analysis::CentreGroup {
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
    pub(crate) fn native(&self) -> molframe::analysis::LeafletOptions {
        molframe::analysis::LeafletOptions {
            connection_distance: self.connection_distance,
            backend: self.backend.into(),
        }
    }

    pub(crate) fn from_native(options: molframe::analysis::LeafletOptions) -> Self {
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
    pub(super) sites: Vec<u32>,
}

pub(crate) fn radial_analysis(
    py: Python<'_>,
    analysis: molframe::Analysis<Vec<molframe::analysis::RadialBin>>,
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
    analysis: molframe::Analysis<Vec<u32>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        Ok(PyList::new(py, values)?.unbind().into_any())
    })
}

pub(crate) fn leaflets_analysis(
    py: Python<'_>,
    analysis: molframe::Analysis<Vec<molframe::analysis::Leaflet>>,
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
            molframe::analysis::radial_distribution(
                structure.positions(),
                &left,
                &right,
                molframe::analysis::RadialDistributionOptions {
                    minimum_distance: options.minimum,
                    maximum_distance: options.maximum,
                    bins: options.bins,
                    volume: options.volume,
                    backend: options.backend.into(),
                },
                periodic_box.as_ref(),
                &crate::core::execution::default_context(),
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
            molframe::analysis::centre_of_mass_radial_distribution(
                structure.positions(),
                &masses,
                &left,
                &right,
                options,
                periodic_box.as_ref(),
                &crate::core::execution::default_context(),
            )
            .map_err(|error| error.to_string())
        })
        .map(|bins| bins.into_iter().map(Into::into).collect())
        .map_err(PyValueError::new_err)
    }

    #[pyo3(signature = (left, right, shell, *, backend=PySpatialBackend::Auto, periodic=false, context=None))]
    fn coordination_numbers(
        &self,
        py: Python<'_>,
        left: &PySelection,
        right: &PySelection,
        shell: (f32, f32),
        backend: PySpatialBackend,
        periodic: bool,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<u32>> {
        let structure = self.structure().clone();
        let left = left.inner.clone();
        let right = right.inner.clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            molframe::analysis::coordination_numbers(
                structure.positions(),
                &left,
                &right,
                molframe::analysis::CoordinationOptions {
                    minimum_distance: shell.0,
                    maximum_distance: shell.1,
                    backend: backend.into(),
                },
                periodic_box.as_ref(),
                &context,
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
            molframe::analysis::identify_leaflets(
                structure.positions(),
                &sites,
                molframe::analysis::LeafletOptions {
                    connection_distance,
                    backend: backend.into(),
                },
                periodic_box.as_ref(),
                &crate::core::execution::default_context(),
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
        molframe::analysis::radial_distribution(
            structure.positions(),
            &left,
            &right,
            options,
            periodic_box.as_ref(),
            &crate::core::execution::default_context(),
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
        molframe::analysis::centre_of_mass_radial_distribution(
            structure.positions(),
            &masses,
            &left,
            &right,
            options,
            periodic_box.as_ref(),
            &crate::core::execution::default_context(),
        )
        .map_err(|error| error.to_string())
    })
    .map(|values| values.into_iter().map(PyRadialBin::from).collect())
    .map_err(PyValueError::new_err)
}

#[pyfunction]
#[pyo3(signature = (structure, left, right, shell, *, backend=PySpatialBackend::Auto, periodic=false, context=None))]
pub(crate) fn coordination_numbers(
    py: Python<'_>,
    structure: &PyStructure,
    left: &PySelection,
    right: &PySelection,
    shell: (f32, f32),
    backend: PySpatialBackend,
    periodic: bool,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<u32>> {
    let structure = structure.structure().clone();
    let left = left.inner.clone();
    let right = right.inner.clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        let periodic_box = periodic_box(&structure, periodic)?;
        molframe::analysis::coordination_numbers(
            structure.positions(),
            &left,
            &right,
            molframe::analysis::CoordinationOptions {
                minimum_distance: shell.0,
                maximum_distance: shell.1,
                backend: backend.into(),
            },
            periodic_box.as_ref(),
            &context,
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
        molframe::analysis::identify_leaflets(
            structure.positions(),
            &sites,
            molframe::analysis::LeafletOptions {
                connection_distance: options.connection_distance,
                backend: options.backend.into(),
            },
            periodic_box.as_ref(),
            &crate::core::execution::default_context(),
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
            let kernel = molframe::analysis::centre_of_mass_radial_distribution_kernel(
                &masses, &left, &right, options,
            );
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
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

fn groups(values: Vec<PySelection>) -> Vec<molframe::analysis::CentreGroup> {
    values
        .into_iter()
        .map(|value| molframe::analysis::CentreGroup { atoms: value.inner })
        .collect()
}

fn rust_options(options: PyRadialOptions) -> molframe::analysis::RadialDistributionOptions {
    molframe::analysis::RadialDistributionOptions {
        minimum_distance: options.minimum,
        maximum_distance: options.maximum,
        bins: options.bins,
        volume: options.volume,
        backend: options.backend.into(),
    }
}
