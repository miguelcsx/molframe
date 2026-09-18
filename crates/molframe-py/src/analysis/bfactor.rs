//! B-factor and TLS validation delegated to native kernels.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::{PyAnalysisPolicy, PySelection};
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "BFactorOutlier", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBFactorOutlier {
    #[pyo3(get)]
    atom: u32,
    #[pyo3(get)]
    value: f64,
    #[pyo3(get)]
    z_score: f64,
}

#[pyclass(name = "BFactorDistribution", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBFactorDistribution {
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    mean: f64,
    #[pyo3(get)]
    variance: f64,
    #[pyo3(get)]
    standard_deviation: f64,
    #[pyo3(get)]
    minimum: f64,
    #[pyo3(get)]
    median: f64,
    #[pyo3(get)]
    maximum: f64,
    #[pyo3(get)]
    outliers: Vec<PyBFactorOutlier>,
}

#[pyclass(name = "TlsModel", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTlsModel(molframe::validate::TlsModel);

#[pymethods]
impl PyTlsModel {
    #[new]
    fn new(
        origin: [f64; 3],
        translation: [[f64; 3]; 3],
        libration: [[f64; 3]; 3],
        screw: [[f64; 3]; 3],
    ) -> Self {
        Self(molframe::validate::TlsModel {
            origin,
            translation,
            libration,
            screw,
        })
    }
}

#[pyclass(name = "TlsGroup", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTlsGroup {
    id: String,
    atoms: PySelection,
    model: PyTlsModel,
}

#[pymethods]
impl PyTlsGroup {
    #[new]
    fn new(id: String, atoms: PySelection, model: PyTlsModel) -> Self {
        Self { id, atoms, model }
    }
}

#[pyclass(name = "TlsBFactorFlag", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTlsBFactorFlag {
    #[pyo3(get)]
    group: String,
    #[pyo3(get)]
    atom: u32,
    #[pyo3(get)]
    observed: f64,
    #[pyo3(get)]
    predicted: f64,
    #[pyo3(get)]
    deviation: f64,
}

#[pyclass(name = "TlsBFactorReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTlsBFactorReport {
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    flags: Vec<PyTlsBFactorFlag>,
}

#[pyfunction]
pub(crate) fn b_factor_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    selection: &PySelection,
    outlier_standard_deviations: f64,
) -> PyResult<PyBFactorDistribution> {
    let structure = structure.structure().clone();
    let selection = selection.inner.clone();
    py.detach(move || {
        molframe::validate::b_factor_distribution(
            &structure,
            &selection,
            outlier_standard_deviations,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn tls_b_factor_consistency(
    py: Python<'_>,
    structure: &PyStructure,
    groups: Vec<PyTlsGroup>,
    maximum_absolute_deviation: f64,
    symmetry_tolerance: f64,
) -> PyResult<PyTlsBFactorReport> {
    let structure = structure.structure().clone();
    let groups = native_groups(groups);
    py.detach(move || {
        molframe::validate::tls_b_factor_consistency(
            &structure,
            &groups,
            maximum_absolute_deviation,
            symmetry_tolerance,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_b_factor_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    selection: &PySelection,
    outlier_standard_deviations: f64,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let selection = selection.inner.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel = molframe::validate::b_factor_distribution_kernel(
                &selection,
                outlier_standard_deviations,
            );
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyBFactorDistribution::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn analyse_tls_b_factor_consistency(
    py: Python<'_>,
    structure: &PyStructure,
    groups: Vec<PyTlsGroup>,
    maximum_absolute_deviation: f64,
    symmetry_tolerance: f64,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let groups = native_groups(groups);
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel = molframe::validate::tls_b_factor_consistency_kernel(
                &groups,
                maximum_absolute_deviation,
                symmetry_tolerance,
            );
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyTlsBFactorReport::from(value)).map(Py::into_any)
    })
}

fn native_groups(groups: Vec<PyTlsGroup>) -> Vec<molframe::validate::TlsGroup> {
    groups
        .into_iter()
        .map(|group| molframe::validate::TlsGroup {
            id: group.id,
            atoms: group.atoms.inner,
            model: group.model.0,
        })
        .collect()
}

impl From<molframe::validate::BFactorOutlier> for PyBFactorOutlier {
    fn from(value: molframe::validate::BFactorOutlier) -> Self {
        Self {
            atom: value.atom.get(),
            value: value.value,
            z_score: value.z_score,
        }
    }
}

impl From<molframe::validate::BFactorDistribution> for PyBFactorDistribution {
    fn from(value: molframe::validate::BFactorDistribution) -> Self {
        Self {
            intended: value.intended,
            assessed: value.assessed,
            mean: value.mean,
            variance: value.variance,
            standard_deviation: value.standard_deviation,
            minimum: value.minimum,
            median: value.median,
            maximum: value.maximum,
            outliers: value.outliers.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<molframe::validate::TlsBFactorReport> for PyTlsBFactorReport {
    fn from(value: molframe::validate::TlsBFactorReport) -> Self {
        Self {
            intended: value.intended,
            assessed: value.assessed,
            flags: value
                .flags
                .into_iter()
                .map(|flag| PyTlsBFactorFlag {
                    group: flag.group,
                    atom: flag.atom.get(),
                    observed: flag.observed,
                    predicted: flag.predicted,
                    deviation: flag.deviation,
                })
                .collect(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
