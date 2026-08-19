//! Chemistry-aware and contact-map comparison delegated to native kernels.

use crate::chemistry::PyComponentDictionary;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::errors::{cad_construction_error, cad_error};
use crate::geometry::borrowed_coordinates;
use crate::query::PyAnalysisPolicy;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "ContactArea", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyContactArea {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
    #[pyo3(get)]
    area: f64,
}

#[pymethods]
impl PyContactArea {
    #[new]
    const fn new(first: u32, second: u32, area: f64) -> Self {
        Self {
            first,
            second,
            area,
        }
    }
}

#[pyclass(name = "CadContact", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCadContact {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
    #[pyo3(get)]
    reference_area: f64,
    #[pyo3(get)]
    model_area: f64,
    #[pyo3(get)]
    lost_area: f64,
}

#[pyclass(name = "LocalCad", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyLocalCad {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    reference_area: f64,
    #[pyo3(get)]
    lost_area: f64,
    #[pyo3(get)]
    score: f64,
}

#[pyclass(name = "CadScore", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCadScore {
    #[pyo3(get)]
    score: f64,
    #[pyo3(get)]
    reference_area: f64,
    #[pyo3(get)]
    lost_area: f64,
    #[pyo3(get)]
    contacts: Vec<PyCadContact>,
    #[pyo3(get)]
    local: Vec<PyLocalCad>,
}

#[pyclass(name = "ContactSimilarity", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyContactSimilarity {
    #[pyo3(get)]
    shared: usize,
    #[pyo3(get)]
    union: usize,
    #[pyo3(get)]
    jaccard: f64,
}

#[pyclass(name = "EquivalentAtomMapping", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEquivalentAtomMapping {
    #[pyo3(get)]
    reference_to_model: Vec<u32>,
}

#[pyclass(name = "LigandRmsd", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLigandRmsd {
    #[pyo3(get)]
    rmsd: f64,
    #[pyo3(get)]
    mapping: PyEquivalentAtomMapping,
}

#[pyfunction]
pub(crate) fn cad_score(
    py: Python<'_>,
    reference: Vec<PyContactArea>,
    model: Vec<PyContactArea>,
) -> PyResult<PyCadScore> {
    let reference = native_areas(reference);
    let model = native_areas(model);
    py.detach(move || pdbiox::compare::cad_score(&reference, &model))
        .map(Into::into)
        .map_err(|error| cad_error(&error))
}

#[pyfunction]
pub(crate) fn cad_contact_areas(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    residues: PyReadonlyArray1<'_, u32>,
    probe: f32,
    density: f32,
) -> PyResult<Vec<PyContactArea>> {
    let positions = borrowed_coordinates(&positions)?;
    let radius_values = radii.as_slice()?;
    let residue_values = residues.as_slice()?;
    py.detach(move || {
        pdbiox::compare::cad_contact_areas(positions, radius_values, residue_values, probe, density)
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(|error| cad_construction_error(&error))
}

#[pyfunction]
pub(crate) fn contact_map_similarity(
    py: Python<'_>,
    first: Vec<(u32, u32)>,
    second: Vec<(u32, u32)>,
) -> PyContactSimilarity {
    py.detach(move || pdbiox::compare::contact_map_similarity(&first, &second))
        .into()
}

#[pyfunction]
pub(crate) fn equivalent_atom_mappings(
    py: Python<'_>,
    dictionary: &PyComponentDictionary,
    component_id: &str,
    limit: usize,
) -> PyResult<Vec<PyEquivalentAtomMapping>> {
    let component = component(dictionary, component_id)?;
    py.detach(move || pdbiox::compare::equivalent_atom_mappings(&component, limit))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn ligand_symmetry_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    dictionary: &PyComponentDictionary,
    component_id: &str,
    limit: usize,
) -> PyResult<PyLigandRmsd> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    let component = component(dictionary, component_id)?;
    py.detach(move || pdbiox::compare::ligand_symmetry_rmsd(reference, model, &component, limit))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_cad_score(
    py: Python<'_>,
    reference: Vec<PyContactArea>,
    model: Vec<PyContactArea>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = native_areas(reference);
    let model = native_areas(model);
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || pdbiox::compare::governed_cad_score(&reference, &model, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyCadScore::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn analyse_cad_contact_areas(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    residues: PyReadonlyArray1<'_, u32>,
    probe: f32,
    density: f32,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let positions = borrowed_coordinates(&positions)?;
    let radius_values = radii.as_slice()?;
    let residue_values = residues.as_slice()?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            pdbiox::compare::governed_cad_contact_areas(
                positions,
                radius_values,
                residue_values,
                probe,
                density,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, values| {
        let items = values
            .into_iter()
            .map(|value| Py::new(py, PyContactArea::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, items)?.unbind().into_any())
    })
}

#[pyfunction]
pub(crate) fn analyse_contact_map_similarity(
    py: Python<'_>,
    first: Vec<(u32, u32)>,
    second: Vec<(u32, u32)>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || pdbiox::compare::governed_contact_map_similarity(&first, &second, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyContactSimilarity::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn analyse_equivalent_atom_mappings(
    py: Python<'_>,
    dictionary: &PyComponentDictionary,
    component_id: &str,
    limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let component = component(dictionary, component_id)?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            pdbiox::compare::governed_equivalent_atom_mappings(&component, limit, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, mapping_list)
}

#[pyfunction]
pub(crate) fn analyse_ligand_symmetry_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    dictionary: &PyComponentDictionary,
    component_id: &str,
    limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    let component = component(dictionary, component_id)?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            pdbiox::compare::governed_ligand_symmetry_rmsd(
                reference, model, &component, limit, &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyLigandRmsd::from(value)).map(Py::into_any)
    })
}

fn component(
    dictionary: &PyComponentDictionary,
    component_id: &str,
) -> PyResult<std::sync::Arc<pdbiox::Component>> {
    pdbiox::ComponentProvider::get(dictionary.0.as_ref(), component_id)
        .map_err(value_error)?
        .ok_or_else(|| PyKeyError::new_err(component_id.to_owned()))
}

fn native_areas(values: Vec<PyContactArea>) -> Vec<pdbiox::compare::ContactArea> {
    values
        .into_iter()
        .map(|value| pdbiox::compare::ContactArea {
            first: value.first,
            second: value.second,
            area: value.area,
        })
        .collect()
}

impl From<pdbiox::compare::ContactArea> for PyContactArea {
    fn from(value: pdbiox::compare::ContactArea) -> Self {
        Self::new(value.first, value.second, value.area)
    }
}

fn mapping_list(
    py: Python<'_>,
    values: Vec<pdbiox::compare::EquivalentAtomMapping>,
) -> PyResult<Py<PyAny>> {
    let items = values
        .into_iter()
        .map(|value| Py::new(py, PyEquivalentAtomMapping::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, items)?.unbind().into_any())
}

impl From<pdbiox::compare::CadContact> for PyCadContact {
    fn from(value: pdbiox::compare::CadContact) -> Self {
        Self {
            first: value.first,
            second: value.second,
            reference_area: value.reference_area,
            model_area: value.model_area,
            lost_area: value.lost_area,
        }
    }
}

impl From<pdbiox::compare::LocalCad> for PyLocalCad {
    fn from(value: pdbiox::compare::LocalCad) -> Self {
        Self {
            residue: value.residue,
            reference_area: value.reference_area,
            lost_area: value.lost_area,
            score: value.score,
        }
    }
}

impl From<pdbiox::compare::CadScore> for PyCadScore {
    fn from(value: pdbiox::compare::CadScore) -> Self {
        Self {
            score: value.score,
            reference_area: value.reference_area,
            lost_area: value.lost_area,
            contacts: value.contacts.into_iter().map(Into::into).collect(),
            local: value.local.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<pdbiox::compare::ContactSimilarity> for PyContactSimilarity {
    fn from(value: pdbiox::compare::ContactSimilarity) -> Self {
        Self {
            shared: value.shared,
            union: value.union,
            jaccard: value.jaccard,
        }
    }
}

impl From<pdbiox::compare::EquivalentAtomMapping> for PyEquivalentAtomMapping {
    fn from(value: pdbiox::compare::EquivalentAtomMapping) -> Self {
        Self {
            reference_to_model: value.reference_to_model.into_vec(),
        }
    }
}

impl From<pdbiox::compare::LigandRmsd> for PyLigandRmsd {
    fn from(value: pdbiox::compare::LigandRmsd) -> Self {
        Self {
            rmsd: value.rmsd,
            mapping: value.mapping.into(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
