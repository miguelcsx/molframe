//! Governed validation reports backed entirely by native kernels.

use super::PyPlanarityOptions;
use crate::chemistry::PyComponentDictionary;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::{PyAnalysisPolicy, PyNamespace, PySelection};
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "AltlocOccupancyOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAltlocOccupancyOptions(pub(crate) molframe::validate::AltlocOccupancyOptions);

#[pymethods]
impl PyAltlocOccupancyOptions {
    #[new]
    fn new(expected_sum: f64, tolerance: f64) -> Self {
        Self(molframe::validate::AltlocOccupancyOptions {
            expected_sum,
            tolerance,
        })
    }
    #[getter]
    fn expected_sum(&self) -> f64 {
        self.0.expected_sum
    }
    #[getter]
    fn tolerance(&self) -> f64 {
        self.0.tolerance
    }
}

#[pyclass(name = "AltlocOccupancyIssue", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAltlocOccupancyIssue {
    MissingOccupancy,
    SumMismatch,
}

#[pyclass(name = "AltlocOccupancyRecord", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAltlocOccupancyRecord {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    atom_name: String,
    #[pyo3(get)]
    alternatives: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    occupancy_sum: Option<f64>,
    #[pyo3(get)]
    issue: Option<PyAltlocOccupancyIssue>,
}

#[pyclass(name = "AltlocOccupancyReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAltlocOccupancyReport {
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    records: Vec<PyAltlocOccupancyRecord>,
}

#[pyclass(name = "ResidueAtomCompleteness", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyResidueAtomCompleteness {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    component: String,
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    missing: Vec<String>,
    #[pyo3(get)]
    ambiguous: Vec<String>,
}

#[pyclass(name = "CcdCompletenessReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCcdCompletenessReport {
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    ambiguous: usize,
    #[pyo3(get)]
    residues: Vec<PyResidueAtomCompleteness>,
}

#[pyclass(name = "PlaneRestraint", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlaneRestraint {
    pub(crate) id: String,
    pub(crate) atoms: PySelection,
}

#[pymethods]
impl PyPlaneRestraint {
    #[new]
    fn new(id: String, atoms: PySelection) -> Self {
        Self { id, atoms }
    }
    #[getter]
    fn id(&self) -> &str {
        &self.id
    }
    #[getter]
    fn atoms(&self) -> PySelection {
        self.atoms.clone()
    }
}

#[pyclass(name = "PlaneRestraintFlag", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlaneRestraintFlag {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    deviation: f64,
}

#[pyclass(name = "PlaneRestraintReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlaneRestraintReport {
    #[pyo3(get)]
    intended: usize,
    #[pyo3(get)]
    assessed: usize,
    #[pyo3(get)]
    flags: Vec<PyPlaneRestraintFlag>,
}

#[pyfunction]
pub(crate) fn validate_altloc_occupancy(
    py: Python<'_>,
    structure: &PyStructure,
    namespace: PyNamespace,
    options: PyAltlocOccupancyOptions,
) -> PyResult<PyAltlocOccupancyReport> {
    let structure = structure.structure().clone();
    py.detach(move || {
        molframe::validate::altloc_occupancy_sums(&structure, namespace.into(), options.0)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_altloc_occupancy(
    py: Python<'_>,
    structure: &PyStructure,
    options: PyAltlocOccupancyOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel = molframe::validate::altloc_occupancy_sums_kernel(options.0);
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, altloc_report_value)
}

#[pyfunction]
pub(crate) fn validate_ccd_completeness(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyCcdCompletenessReport> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    py.detach(move || {
        molframe::validate::ccd_missing_atoms(&structure, dictionary.as_ref(), &policy)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_ccd_completeness(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel = molframe::validate::ccd_missing_atoms_kernel(dictionary.as_ref());
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, ccd_report_value)
}

#[pyfunction]
pub(crate) fn validate_plane_restraints(
    py: Python<'_>,
    structure: &PyStructure,
    restraints: Vec<PyPlaneRestraint>,
    options: PyPlanarityOptions,
) -> PyResult<PyPlaneRestraintReport> {
    let structure = structure.structure().clone();
    let restraints = rust_restraints(restraints);
    py.detach(move || {
        molframe::validate::plane_restraint_outliers(&structure, &restraints, options.0)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_plane_restraints(
    py: Python<'_>,
    structure: &PyStructure,
    restraints: Vec<PyPlaneRestraint>,
    options: PyPlanarityOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let restraints = rust_restraints(restraints);
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            let kernel =
                molframe::validate::plane_restraint_outliers_kernel(&restraints, options.0);
            molframe::analysis::analyse_structure(
                &structure,
                &policy,
                &kernel,
                &crate::core::execution::default_context(),
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, plane_report_value)
}

fn rust_restraints(values: Vec<PyPlaneRestraint>) -> Vec<molframe::validate::PlaneRestraint> {
    values
        .into_iter()
        .map(|value| molframe::validate::PlaneRestraint {
            id: value.id,
            atoms: value.atoms.inner,
        })
        .collect()
}

fn altloc_report_value(
    py: Python<'_>,
    value: molframe::validate::AltlocOccupancyReport,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyAltlocOccupancyReport::from(value))?.into_any())
}
fn ccd_report_value(
    py: Python<'_>,
    value: molframe::validate::CcdCompletenessReport,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyCcdCompletenessReport::from(value))?.into_any())
}
fn plane_report_value(
    py: Python<'_>,
    value: molframe::validate::PlaneRestraintReport,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyPlaneRestraintReport::from(value))?.into_any())
}

impl From<molframe::validate::AltlocOccupancyIssue> for PyAltlocOccupancyIssue {
    fn from(value: molframe::validate::AltlocOccupancyIssue) -> Self {
        match value {
            molframe::validate::AltlocOccupancyIssue::MissingOccupancy => Self::MissingOccupancy,
            molframe::validate::AltlocOccupancyIssue::SumMismatch => Self::SumMismatch,
        }
    }
}
impl From<molframe::validate::AltlocOccupancyRecord> for PyAltlocOccupancyRecord {
    fn from(value: molframe::validate::AltlocOccupancyRecord) -> Self {
        Self {
            residue: value.residue.get(),
            atom_name: value.atom_name,
            alternatives: value.alternatives,
            assessed: value.assessed,
            occupancy_sum: value.occupancy_sum,
            issue: value.issue.map(Into::into),
        }
    }
}
impl From<molframe::validate::AltlocOccupancyReport> for PyAltlocOccupancyReport {
    fn from(value: molframe::validate::AltlocOccupancyReport) -> Self {
        Self {
            intended: value.intended,
            assessed: value.assessed,
            records: value.records.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<molframe::validate::ResidueAtomCompleteness> for PyResidueAtomCompleteness {
    fn from(value: molframe::validate::ResidueAtomCompleteness) -> Self {
        Self {
            residue: value.residue.get(),
            component: value.component,
            intended: value.intended,
            assessed: value.assessed,
            missing: value.missing,
            ambiguous: value.ambiguous,
        }
    }
}
impl From<molframe::validate::CcdCompletenessReport> for PyCcdCompletenessReport {
    fn from(value: molframe::validate::CcdCompletenessReport) -> Self {
        Self {
            intended: value.intended,
            assessed: value.assessed,
            ambiguous: value.ambiguous,
            residues: value.residues.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<molframe::validate::PlaneRestraintFlag> for PyPlaneRestraintFlag {
    fn from(value: molframe::validate::PlaneRestraintFlag) -> Self {
        Self {
            id: value.id,
            deviation: value.deviation,
        }
    }
}
impl From<molframe::validate::PlaneRestraintReport> for PyPlaneRestraintReport {
    fn from(value: molframe::validate::PlaneRestraintReport) -> Self {
        Self {
            intended: value.intended,
            assessed: value.assessed,
            flags: value.flags.into_iter().map(Into::into).collect(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
