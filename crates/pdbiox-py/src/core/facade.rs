//! Direct bindings for public facade projections that are not structure methods.

use crate::chemistry::PyComponentDictionary;
use crate::geometry::PyBackboneTorsions;
use crate::graph::PySpatialBackend;
use crate::io::PyLimits;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "BondInference", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBondInference(pub(crate) pdbiox::BondInference);

#[pymethods]
impl PyBondInference {
    #[new]
    #[pyo3(signature = (scale=1.15, lower_bound=0.4, *, exclude_across_chains=false, respect_existing=true, backend=PySpatialBackend::Auto))]
    fn new(
        scale: f32,
        lower_bound: f32,
        exclude_across_chains: bool,
        respect_existing: bool,
        backend: PySpatialBackend,
    ) -> Self {
        Self(pdbiox::BondInference {
            scale,
            lower_bound,
            exclude_across_chains,
            respect_existing,
            backend: backend.into(),
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::BondInference::default())
    }

    #[getter]
    const fn scale(&self) -> f32 {
        self.0.scale
    }

    #[getter]
    const fn lower_bound(&self) -> f32 {
        self.0.lower_bound
    }

    #[getter]
    const fn exclude_across_chains(&self) -> bool {
        self.0.exclude_across_chains
    }

    #[getter]
    const fn respect_existing(&self) -> bool {
        self.0.respect_existing
    }

    #[getter]
    fn backend(&self) -> PySpatialBackend {
        self.0.backend.into()
    }
}

#[pyclass(name = "BondInferenceReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondInferenceReport {
    structure: PyStructure,
    #[pyo3(get)]
    skipped_atoms: Vec<u32>,
}

impl PyBondInferenceReport {
    pub(crate) fn from_native(report: pdbiox::BondInferenceReport) -> Self {
        Self {
            structure: PyStructure::new(report.structure),
            skipped_atoms: report
                .skipped_atoms
                .into_iter()
                .map(pdbiox::AtomIndex::get)
                .collect(),
        }
    }
}

#[pymethods]
impl PyBondInferenceReport {
    #[getter]
    fn structure(&self) -> PyStructure {
        self.structure.clone()
    }
}

#[pyclass(name = "ProteinAlphaTrace", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyProteinAlphaTrace {
    #[pyo3(get)]
    chain: u32,
    #[pyo3(get)]
    positions: Vec<Option<[f32; 3]>>,
}

#[pyclass(name = "BackboneTorsionRecord", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBackboneTorsionRecord {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    torsions: PyBackboneTorsions,
}

#[pyclass(name = "SideChainTorsionRecord", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySideChainTorsionRecord {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    atoms: Vec<String>,
    #[pyo3(get)]
    torsions: Vec<Option<f64>>,
}

#[pyclass(name = "SideChainTorsionReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySideChainTorsionReport {
    #[pyo3(get)]
    records: Vec<PySideChainTorsionRecord>,
    #[pyo3(get)]
    findings: Vec<String>,
    #[pyo3(get)]
    dictionary_version: String,
}

#[pyfunction]
pub(crate) fn default_limits(py: Python<'_>) -> PyLimits {
    py.detach(move || -> PyLimits { PyLimits::from_inner(pdbiox::default_limits()) })
}

#[pyfunction]
pub(crate) fn infer_bonds(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyBondInference,
) -> PyResult<PyBondInferenceReport> {
    let structure = structure.structure().clone();
    let options = options.0;
    py.detach(move || {
        pdbiox::infer_bonds(
            &structure,
            options,
            &crate::core::execution::default_context(),
        )
    })
    .map(PyBondInferenceReport::from_native)
    .map_err(|finding| crate::errors::read_error(py, std::slice::from_ref(&finding)))
}

#[pyfunction]
pub(crate) fn structure_protein_alpha_traces(
    py: Python<'_>,
    structure: &PyStructure,
) -> PyResult<Vec<PyProteinAlphaTrace>> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::structure_protein_alpha_traces(&structure))
        .map(|traces| {
            traces
                .into_iter()
                .map(|trace| PyProteinAlphaTrace {
                    chain: trace.chain.get(),
                    positions: trace.positions,
                })
                .collect()
        })
        .map_err(|finding| crate::errors::read_error(py, std::slice::from_ref(&finding)))
}

#[pyfunction]
pub(crate) fn structure_backbone_torsions(
    py: Python<'_>,
    structure: &PyStructure,
) -> PyResult<Vec<PyBackboneTorsionRecord>> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::structure_backbone_torsions(&structure))
        .map(|records| {
            records
                .into_iter()
                .map(PyBackboneTorsionRecord::from)
                .collect()
        })
        .map_err(|finding| crate::errors::read_error(py, std::slice::from_ref(&finding)))
}

#[pyfunction]
pub(crate) fn structure_backbone_torsions_model(
    py: Python<'_>,
    structure: &PyStructure,
    model: u32,
) -> PyResult<Vec<PyBackboneTorsionRecord>> {
    let structure = structure.structure().clone();
    py.detach(move || {
        pdbiox::structure_backbone_torsions_model(&structure, pdbiox::ModelIndex::new(model))
    })
    .map(|records| {
        records
            .into_iter()
            .map(PyBackboneTorsionRecord::from)
            .collect()
    })
    .map_err(|finding| crate::errors::read_error(py, std::slice::from_ref(&finding)))
}

#[pyfunction]
pub(crate) fn structure_side_chain_torsions(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    policy: &PyAnalysisPolicy,
) -> PyResult<PySideChainTorsionReport> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    py.detach(move || {
        pdbiox::structure_side_chain_torsions(&structure, dictionary.as_ref(), &policy)
    })
    .map(|report| PySideChainTorsionReport {
        records: report
            .records
            .into_iter()
            .map(|record| PySideChainTorsionRecord {
                residue: record.residue.get(),
                atoms: record.atoms.into_iter().map(String::from).collect(),
                torsions: record.torsions.into_vec(),
            })
            .collect(),
        findings: report
            .findings
            .into_iter()
            .map(|finding| finding.to_string())
            .collect(),
        dictionary_version: report.dictionary_version.as_str().to_owned(),
    })
    .map_err(|finding| crate::errors::read_error(py, std::slice::from_ref(&finding)))
}

impl From<pdbiox::BackboneTorsionRecord> for PyBackboneTorsionRecord {
    fn from(value: pdbiox::BackboneTorsionRecord) -> Self {
        Self {
            residue: value.residue.get(),
            torsions: value.torsions.into(),
        }
    }
}

#[cfg(test)]
#[path = "facade_tests.rs"]
mod tests;
