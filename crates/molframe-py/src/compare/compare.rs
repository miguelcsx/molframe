//! Vectorized structure scores delegated to `molframe-compare`.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::core::execution::PyExecutionContext;
use crate::errors::ce_error;
use crate::geometry::borrowed_coordinates;
use crate::query::{PyAnalysisPolicy, PyNamespace};
use crate::structure::PyStructure;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyFloat;

macro_rules! score {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name(
            py: Python<'_>,
            model: PyReadonlyArray2<'_, f32>,
            reference: PyReadonlyArray2<'_, f32>,
        ) -> PyResult<f64> {
            let model = borrowed_coordinates(&model)?;
            let reference = borrowed_coordinates(&reference)?;
            py.detach(|| $kernel(model, reference)).map_err(value_error)
        }
    };
}

score!(tm_score, molframe::compare::tm_score);
score!(gdt_ts, molframe::compare::gdt_ts);
score!(gdt_ha, molframe::compare::gdt_ha);

#[pyclass(name = "EmptyLddtPolicy", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyEmptyLddtPolicy {
    Perfect,
    Error,
}

#[pyclass(name = "LddtOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLddtOptions(molframe::compare::LddtOptions);

#[pymethods]
impl PyLddtOptions {
    #[new]
    #[pyo3(signature = (
        inclusion_radius,
        minimum_reference_distance,
        tolerances,
        empty_policy,
    ))]
    fn new(
        inclusion_radius: f64,
        minimum_reference_distance: f64,
        tolerances: Vec<f64>,
        empty_policy: PyEmptyLddtPolicy,
    ) -> Self {
        Self(molframe::compare::LddtOptions {
            inclusion_radius,
            minimum_reference_distance,
            tolerances: tolerances.into_boxed_slice(),
            empty_policy: empty_policy.into(),
        })
    }

    #[staticmethod]
    fn standard(inclusion_radius: f64) -> Self {
        Self(molframe::compare::LddtOptions::standard(inclusion_radius))
    }
}

#[pyfunction(signature = (model, reference, options, *, context=None))]
pub(crate) fn lddt(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    options: &PyLddtOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<f64> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    let options = options.0.clone();
    let execution = context.map_or_else(molframe::core::ExecutionContext::default, |value| {
        value.native()
    });
    py.detach(|| molframe::compare::lddt_with_options(model, reference, &options, &execution))
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn gdt(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    cutoffs: Vec<f64>,
) -> PyResult<f64> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    let cutoffs = cutoffs.into_boxed_slice();
    py.detach(|| molframe::compare::gdt_with_cutoffs(model, reference, &cutoffs))
        .map_err(value_error)
}

#[pyclass(name = "CeSignificanceProfile", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCeSignificanceProfile {
    OriginalWindowEight,
}

#[pyclass(name = "CeOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCeOptions(molframe::compare::CeOptions);

impl PyCeOptions {
    pub(crate) const fn inner(self) -> molframe::compare::CeOptions {
        self.0
    }
}

#[pymethods]
impl PyCeOptions {
    #[new]
    #[pyo3(signature = (
        window_size,
        max_gap,
        max_paths,
        fragment_similarity_threshold,
        path_similarity_threshold,
        significance,
        memory_limit_bytes=100_000_000
    ))]
    fn new(
        window_size: usize,
        max_gap: usize,
        max_paths: usize,
        fragment_similarity_threshold: f64,
        path_similarity_threshold: f64,
        significance: Option<PyCeSignificanceProfile>,
        memory_limit_bytes: usize,
    ) -> Self {
        Self(molframe::compare::CeOptions {
            window_size,
            max_gap,
            max_paths,
            fragment_similarity_threshold,
            path_similarity_threshold,
            significance: significance.map(Into::into),
            memory_limit_bytes,
        })
    }

    #[staticmethod]
    fn original() -> Self {
        Self(molframe::compare::CeOptions::original())
    }
}

#[pyclass(name = "CeAlignment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCeAlignment {
    #[pyo3(get)]
    reference_indices: Vec<usize>,
    #[pyo3(get)]
    mobile_indices: Vec<usize>,
    #[pyo3(get)]
    fragment_count: usize,
    #[pyo3(get)]
    similarity: f64,
    #[pyo3(get)]
    z_score: Option<f64>,
    #[pyo3(get)]
    rmsd: f64,
}

#[pyfunction]
pub(crate) fn ce_align(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    mobile: PyReadonlyArray2<'_, f32>,
    options: PyCeOptions,
) -> PyResult<PyCeAlignment> {
    let reference = borrowed_coordinates(&reference)?;
    let mobile = borrowed_coordinates(&mobile)?;
    py.detach(|| molframe::compare::ce_align(reference, mobile, options.0))
        .map(Into::into)
        .map_err(|error| ce_error(&error))
}

#[pyfunction]
pub(crate) fn ce_alignments(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    mobile: PyReadonlyArray2<'_, f32>,
    options: PyCeOptions,
) -> PyResult<Vec<PyCeAlignment>> {
    let reference = borrowed_coordinates(&reference)?;
    let mobile = borrowed_coordinates(&mobile)?;
    py.detach(|| molframe::compare::ce_alignments(reference, mobile, options.0))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(|error| ce_error(&error))
}

#[pyclass(name = "EmptyQsPolicy", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyEmptyQsPolicy {
    Perfect,
    Error,
}

#[pyclass(name = "QsOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyQsOptions(molframe::compare::QsOptions);

#[pymethods]
impl PyQsOptions {
    #[new]
    fn new(contact_distance: f32, empty_policy: PyEmptyQsPolicy) -> Self {
        Self(molframe::compare::QsOptions {
            contact_distance,
            empty_policy: empty_policy.into(),
        })
    }

    #[staticmethod]
    fn standard(contact_distance: f32) -> Self {
        Self(molframe::compare::QsOptions::standard(contact_distance))
    }
}

#[pyclass(name = "DockQOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDockQOptions(molframe::compare::DockQOptions);

#[pymethods]
impl PyDockQOptions {
    #[new]
    fn new(contact_distance: f32, ligand_scale: f64, interface_scale: f64) -> Self {
        Self(molframe::compare::DockQOptions {
            contact_distance,
            ligand_scale,
            interface_scale,
        })
    }
}

#[pyclass(name = "DockQ", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDockQ {
    #[pyo3(get)]
    fnat: f64,
    #[pyo3(get)]
    ligand_rmsd: f64,
    #[pyo3(get)]
    interface_rmsd: f64,
    #[pyo3(get)]
    score: f64,
}

#[pyfunction]
pub(crate) fn qs_score(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    namespace: PyNamespace,
    options: PyQsOptions,
) -> PyResult<f64> {
    let model = model.structure().clone();
    let native = native.structure().clone();
    let first_chain = first_chain.to_owned();
    let second_chain = second_chain.to_owned();
    py.detach(|| {
        molframe::compare::qs_score_in_namespace(
            &model,
            &native,
            &first_chain,
            &second_chain,
            namespace.into(),
            options.0,
        )
    })
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn dockq(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    receptor: &str,
    ligand: &str,
    namespace: PyNamespace,
    options: PyDockQOptions,
) -> PyResult<PyDockQ> {
    let model = model.structure().clone();
    let native = native.structure().clone();
    let receptor = receptor.to_owned();
    let ligand = ligand.to_owned();
    py.detach(|| {
        molframe::compare::dockq_in_namespace(
            &model,
            &native,
            &receptor,
            &ligand,
            namespace.into(),
            options.0,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn weighted_rmsd(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    weights: Vec<f64>,
) -> PyResult<f64> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    py.detach(move || molframe::compare::weighted_rmsd(model, reference, &weights))
        .map_err(value_error)
}

#[pyfunction(signature = (model, reference, options, policy, *, context=None))]
pub(crate) fn analyse_lddt(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    options: &PyLddtOptions,
    policy: &PyAnalysisPolicy,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    let options = options.0.clone();
    let policy = policy.inner.clone();
    let execution = context.map_or_else(molframe::core::ExecutionContext::default, |value| {
        value.native()
    });
    let analysis = py
        .detach(|| {
            molframe::compare::governed_lddt(model, reference, &options, &policy, &execution)
        })
        .map_err(value_error)?;
    scalar_analysis(py, analysis)
}

macro_rules! governed_score {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name(
            py: Python<'_>,
            model: PyReadonlyArray2<'_, f32>,
            reference: PyReadonlyArray2<'_, f32>,
            policy: &PyAnalysisPolicy,
        ) -> PyResult<PyAnalysis> {
            let model = borrowed_coordinates(&model)?;
            let reference = borrowed_coordinates(&reference)?;
            let policy = policy.inner.clone();
            let analysis = py
                .detach(|| $kernel(model, reference, &policy))
                .map_err(value_error)?;
            scalar_analysis(py, analysis)
        }
    };
}

governed_score!(analyse_tm_score, molframe::compare::governed_tm_score);
governed_score!(analyse_gdt_ts, molframe::compare::governed_gdt_ts);
governed_score!(analyse_gdt_ha, molframe::compare::governed_gdt_ha);

#[pyfunction]
pub(crate) fn analyse_dockq(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    receptor: &str,
    ligand: &str,
    options: PyDockQOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let model = model.structure().clone();
    let native = native.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            molframe::compare::governed_dockq(&model, &native, receptor, ligand, options.0, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyDockQ::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn analyse_qs_score(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    options: PyQsOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let model = model.structure().clone();
    let native = native.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            molframe::compare::governed_qs_score(
                &model,
                &native,
                first_chain,
                second_chain,
                options.0,
                &policy,
            )
        })
        .map_err(value_error)?;
    scalar_analysis(py, analysis)
}

fn scalar_analysis(py: Python<'_>, analysis: molframe::Analysis<f64>) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        Ok(PyFloat::new(py, value).unbind().into_any())
    })
}

impl From<PyEmptyLddtPolicy> for molframe::compare::EmptyLddtPolicy {
    fn from(value: PyEmptyLddtPolicy) -> Self {
        match value {
            PyEmptyLddtPolicy::Perfect => Self::Perfect,
            PyEmptyLddtPolicy::Error => Self::Error,
        }
    }
}

impl From<PyCeSignificanceProfile> for molframe::compare::CeSignificanceProfile {
    fn from(value: PyCeSignificanceProfile) -> Self {
        match value {
            PyCeSignificanceProfile::OriginalWindowEight => Self::OriginalWindowEight,
        }
    }
}

impl From<PyEmptyQsPolicy> for molframe::compare::EmptyQsPolicy {
    fn from(value: PyEmptyQsPolicy) -> Self {
        match value {
            PyEmptyQsPolicy::Perfect => Self::Perfect,
            PyEmptyQsPolicy::Error => Self::Error,
        }
    }
}

impl From<molframe::compare::CeAlignment> for PyCeAlignment {
    fn from(value: molframe::compare::CeAlignment) -> Self {
        Self {
            reference_indices: value.reference_indices,
            mobile_indices: value.mobile_indices,
            fragment_count: value.fragment_count,
            similarity: value.similarity,
            z_score: value.z_score,
            rmsd: value.rmsd,
        }
    }
}

impl From<molframe::compare::DockQ> for PyDockQ {
    fn from(value: molframe::compare::DockQ) -> Self {
        Self {
            fnat: value.fnat,
            ligand_rmsd: value.ligand_rmsd,
            interface_rmsd: value.interface_rmsd,
            score: value.score,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
