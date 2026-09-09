//! Fragment-library and elastic-network bindings.

use crate::core::parallel::PyReductionPolicy;
use crate::geometry::coordinates;
use crate::graph::PySpatialBackend;
use crate::query::PySelection;
use crate::structure::PyStructure;
use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::periodic::periodic_box;

#[pyclass(name = "FragmentReference", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFragmentReference(pdbiox::analysis::FragmentReference);

#[pymethods]
impl PyFragmentReference {
    #[new]
    fn new(id: &str, coordinates_array: PyReadonlyArray2<'_, f32>) -> PyResult<Self> {
        Ok(Self(pdbiox::analysis::FragmentReference {
            id: id.into(),
            coordinates: coordinates(coordinates_array)?,
        }))
    }
}

#[pyclass(name = "FragmentMatch", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFragmentMatch {
    #[pyo3(get)]
    start: usize,
    #[pyo3(get)]
    fragment_id: String,
    #[pyo3(get)]
    rmsd: f64,
}

#[pyclass(name = "GnmOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyGnmOptions(pdbiox::analysis::GnmOptions);

#[pymethods]
impl PyGnmOptions {
    #[new]
    #[pyo3(signature = (
        contact_distance,
        mode_count,
        zero_mode_tolerance,
        memory_limit_bytes,
        backend,
        reduction=PyReductionPolicy::Deterministic,
    ))]
    fn new(
        contact_distance: f32,
        mode_count: usize,
        zero_mode_tolerance: f64,
        memory_limit_bytes: usize,
        backend: PySpatialBackend,
        reduction: PyReductionPolicy,
    ) -> Self {
        Self(pdbiox::analysis::GnmOptions {
            contact_distance,
            mode_count,
            zero_mode_tolerance,
            memory_limit_bytes,
            backend: backend.into(),
            reduction: reduction.into(),
        })
    }

    /// Whether the eigensolver may reorder its floating-point reductions.
    #[getter]
    fn reduction(&self) -> PyReductionPolicy {
        PyReductionPolicy::from(self.0.reduction)
    }
}

impl PyGnmOptions {
    pub(crate) const fn from_native(value: pdbiox::analysis::GnmOptions) -> Self {
        Self(value)
    }

    pub(crate) const fn native(self) -> pdbiox::analysis::GnmOptions {
        self.0
    }
}

#[pyclass(name = "GaussianNetworkModel", frozen, skip_from_py_object)]
pub(crate) struct PyGaussianNetworkModel {
    sites: Py<PyArray1<u32>>,
    eigenvalues: Py<PyArray1<f64>>,
    modes: Py<PyArray2<f64>>,
    #[pyo3(get)]
    zero_modes: usize,
}

#[pymethods]
impl PyGaussianNetworkModel {
    #[getter]
    fn sites(&self, py: Python<'_>) -> Py<PyArray1<u32>> {
        self.sites.clone_ref(py)
    }

    #[getter]
    fn eigenvalues(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.eigenvalues.clone_ref(py)
    }

    #[getter]
    fn modes(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.modes.clone_ref(py)
    }
}

#[pyfunction]
pub(crate) fn map_fragments(
    py: Python<'_>,
    trace: Vec<Option<[f32; 3]>>,
    library: Vec<PyFragmentReference>,
    maximum_rmsd: f64,
) -> PyResult<Vec<PyFragmentMatch>> {
    let library = library.into_iter().map(|value| value.0).collect::<Vec<_>>();
    py.detach(move || pdbiox::analysis::map_fragments(&trace, &library, maximum_rmsd))
        .map(|values| values.into_iter().map(PyFragmentMatch::from).collect())
        .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, sites, options, *, periodic=false, context=None))]
pub(crate) fn gaussian_network_model(
    py: Python<'_>,
    structure: &PyStructure,
    sites: &PySelection,
    options: PyGnmOptions,
    periodic: bool,
    context: Option<&crate::core::execution::PyExecutionContext>,
) -> PyResult<PyGaussianNetworkModel> {
    let structure = structure.structure().clone();
    let sites = sites.inner.clone();
    let context = context.map_or_else(
        crate::core::execution::default_context,
        crate::core::execution::PyExecutionContext::native,
    );
    let value = py
        .detach(move || {
            let periodic_box = periodic_box(&structure, periodic)?;
            pdbiox::analysis::gaussian_network_model(
                structure.positions(),
                &sites,
                options.0,
                periodic_box.as_ref(),
                &context,
            )
            .map_err(|error| error.to_string())
        })
        .map_err(PyValueError::new_err)?;
    PyGaussianNetworkModel::new(py, value)
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (sites, options, *, periodic=false, context=None))]
    fn gaussian_network_model(
        &self,
        py: Python<'_>,
        sites: &PySelection,
        options: PyGnmOptions,
        periodic: bool,
        context: Option<&crate::core::execution::PyExecutionContext>,
    ) -> PyResult<PyGaussianNetworkModel> {
        let structure = self.structure().clone();
        let sites = sites.inner.clone();
        let context = context.map_or_else(
            crate::core::execution::default_context,
            crate::core::execution::PyExecutionContext::native,
        );
        let value = py
            .detach(move || {
                let periodic_box = periodic_box(&structure, periodic)?;
                pdbiox::analysis::gaussian_network_model(
                    structure.positions(),
                    &sites,
                    options.0,
                    periodic_box.as_ref(),
                    &context,
                )
                .map_err(|error| error.to_string())
            })
            .map_err(PyValueError::new_err)?;
        PyGaussianNetworkModel::new(py, value)
    }
}

impl PyGaussianNetworkModel {
    pub(crate) fn new(
        py: Python<'_>,
        value: pdbiox::analysis::GaussianNetworkModel,
    ) -> PyResult<Self> {
        let rows = value.modes.len();
        let columns = value.sites.len();
        let modes = value.modes.into_iter().flatten().collect::<Vec<_>>();
        let modes = Array2::from_shape_vec((rows, columns), modes).map_err(value_error)?;
        Ok(Self {
            sites: Array1::from_vec(value.sites).into_pyarray(py).unbind(),
            eigenvalues: Array1::from_vec(value.eigenvalues)
                .into_pyarray(py)
                .unbind(),
            modes: modes.into_pyarray(py).unbind(),
            zero_modes: value.zero_modes,
        })
    }
}

impl From<pdbiox::analysis::FragmentMatch> for PyFragmentMatch {
    fn from(value: pdbiox::analysis::FragmentMatch) -> Self {
        Self {
            start: value.start,
            fragment_id: value.fragment_id.into(),
            rmsd: value.rmsd,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
