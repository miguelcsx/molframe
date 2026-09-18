//! Owned projections for CCD components and chemistry-provider results.

use super::{PyComponentDictionary, PyElement, PyPeoeOptions};
use crate::bonds::PyBondOrder;
use crate::contract::{PyCoverage, PyDiagnostic};
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::exceptions::{PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::sync::{Arc, RwLock};

#[pyclass(name = "ComponentKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyComponentKind {
    AminoAcid,
    Nucleotide,
    Saccharide,
    Lipid,
    NonPolymer,
    Solvent,
    Ion,
    Unknown,
}

#[pyclass(name = "StereoConfiguration", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyStereoConfiguration {
    R,
    S,
    Mixed,
}

#[pyclass(name = "ComponentAtom", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentAtom {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    alternate_name: Option<String>,
    #[pyo3(get)]
    element: PyElement,
    #[pyo3(get)]
    charge: i8,
    #[pyo3(get)]
    aromatic: bool,
    #[pyo3(get)]
    leaving: bool,
    #[pyo3(get)]
    stereo: Option<PyStereoConfiguration>,
}

#[pyclass(name = "ComponentBond", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentBond {
    #[pyo3(get)]
    atom_a: String,
    #[pyo3(get)]
    atom_b: String,
    #[pyo3(get)]
    order: PyBondOrder,
    #[pyo3(get)]
    aromatic: bool,
    #[pyo3(get)]
    stereo: Option<PyStereoConfiguration>,
}

#[pyclass(name = "Component", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponent(pub(crate) Arc<molframe::Component>);

#[pymethods]
impl PyComponent {
    #[getter]
    fn id(&self) -> String {
        self.0.id.to_string()
    }
    #[getter]
    fn name(&self) -> String {
        self.0.name.to_string()
    }
    #[getter]
    fn kind(&self) -> PyComponentKind {
        self.0.kind.into()
    }
    #[getter]
    fn parent(&self) -> Option<String> {
        self.0.parent.as_deref().map(str::to_owned)
    }
    #[getter]
    fn one_letter_code(&self) -> Option<String> {
        self.0
            .one_letter_code
            .map(char::from)
            .map(|value| value.to_string())
    }
    #[getter]
    fn formula(&self) -> Option<String> {
        self.0.formula.as_deref().map(str::to_owned)
    }
    #[getter]
    fn atoms(&self) -> Vec<PyComponentAtom> {
        self.0.atoms.iter().cloned().map(Into::into).collect()
    }
    #[getter]
    fn bonds(&self) -> Vec<PyComponentBond> {
        self.0.bonds.iter().cloned().map(Into::into).collect()
    }
    #[getter]
    fn ideal_coordinates(&self) -> Option<Vec<[f32; 3]>> {
        self.0
            .ideal_coordinates
            .as_deref()
            .map(<[[f32; 3]]>::to_vec)
    }
    #[getter]
    fn model_coordinates(&self) -> Option<Vec<[f32; 3]>> {
        self.0
            .model_coordinates
            .as_deref()
            .map(<[[f32; 3]]>::to_vec)
    }
    fn atom(&self, name: &str) -> Option<PyComponentAtom> {
        self.0.atom(name).cloned().map(Into::into)
    }
}

#[pyclass(name = "ComponentCoverage", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentCoverage {
    #[pyo3(get)]
    coverage: PyCoverage,
    #[pyo3(get)]
    findings: Vec<PyDiagnostic>,
    #[pyo3(get)]
    dictionary_version: String,
}

#[pyclass(name = "EquivalenceClasses", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEquivalenceClasses(pub(crate) molframe::EquivalenceClasses);

#[pymethods]
impl PyEquivalenceClasses {
    #[getter]
    fn classes(&self) -> Vec<Vec<u32>> {
        self.0
            .classes()
            .iter()
            .map(|class| class.iter().copied().collect())
            .collect()
    }

    fn class_of(&self, atom: u32) -> Option<u32> {
        self.0.class_of(atom)
    }

    fn equivalent(&self, atom_a: u32, atom_b: u32) -> bool {
        self.0.equivalent(atom_a, atom_b)
    }
}

#[pyclass(name = "EquivalenceCache", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyEquivalenceCache {
    inner: RwLock<molframe::EquivalenceCache>,
}

#[pymethods]
impl PyEquivalenceCache {
    #[new]
    fn new(version: &str) -> Self {
        Self {
            inner: RwLock::new(molframe::EquivalenceCache::new(
                molframe::DictionaryVersion::new(version),
            )),
        }
    }

    #[getter]
    fn version(&self) -> PyResult<String> {
        self.inner
            .read()
            .map(|cache| cache.version().as_str().to_owned())
            .map_err(|_| PyRuntimeError::new_err("equivalence cache lock is poisoned"))
    }

    #[getter]
    fn len(&self) -> PyResult<usize> {
        self.inner
            .read()
            .map(|cache| cache.len())
            .map_err(|_| PyRuntimeError::new_err("equivalence cache lock is poisoned"))
    }

    #[getter]
    fn is_empty(&self) -> PyResult<bool> {
        self.inner
            .read()
            .map(|cache| cache.is_empty())
            .map_err(|_| PyRuntimeError::new_err("equivalence cache lock is poisoned"))
    }

    fn get(&self, component: &PyComponent) -> PyResult<PyEquivalenceClasses> {
        let classes = self
            .inner
            .write()
            .map_err(|_| PyRuntimeError::new_err("equivalence cache lock is poisoned"))?
            .get(&component.0)
            .clone();
        Ok(PyEquivalenceClasses((*classes).clone()))
    }
}

#[pyclass(name = "PolymerAtomRole", frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PyPolymerAtomRole(pub(crate) molframe::PolymerAtomRole);

#[pymethods]
impl PyPolymerAtomRole {
    #[staticmethod]
    fn from_code(code: i64) -> Option<Self> {
        molframe::PolymerAtomRole::from_code(code).map(Self)
    }
    #[getter]
    fn code(&self) -> i64 {
        self.0.code()
    }
    fn intersects(&self, required: &Self) -> bool {
        self.0.intersects(required.0)
    }
    fn union(&self, other: &Self) -> Self {
        Self(self.0.union(other.0))
    }
    #[classattr]
    fn unknown() -> Self {
        Self(molframe::PolymerAtomRole::UNKNOWN)
    }
    #[classattr]
    fn protein_backbone() -> Self {
        Self(molframe::PolymerAtomRole::PROTEIN_BACKBONE)
    }
    #[classattr]
    fn nucleic_backbone() -> Self {
        Self(molframe::PolymerAtomRole::NUCLEIC_BACKBONE)
    }
}

#[pymethods]
impl PyComponentDictionary {
    fn get(&self, component_id: &str) -> PyResult<Option<PyComponent>> {
        molframe::ComponentProvider::get(self.0.as_ref(), component_id)
            .map(|component| component.map(PyComponent))
            .map_err(value_error)
    }

    fn equivalence_classes(&self, component_id: &str) -> PyResult<Option<PyEquivalenceClasses>> {
        let component = self.get(component_id)?;
        Ok(component
            .map(|component| PyEquivalenceClasses(molframe::equivalence_classes(&component.0))))
    }

    pub(crate) fn coverage(
        &self,
        structure: &PyStructure,
        policy: &PyAnalysisPolicy,
    ) -> PyResult<PyComponentCoverage> {
        molframe::component_coverage(structure.structure(), self.0.as_ref(), &policy.inner)
            .map(Into::into)
            .map_err(value_error)
    }

    fn side_chain_definition(
        &self,
        component_id: &str,
        nitrogen: &str,
        alpha_carbon: &str,
        side_chain_atoms: &Bound<'_, pyo3::types::PyList>,
    ) -> PyResult<Option<Vec<String>>> {
        let side_chain_atoms = side_chain_atoms.extract::<Vec<String>>()?;
        let Some(component) = self.get(component_id)? else {
            return Err(PyKeyError::new_err(component_id.to_owned()));
        };
        let atom_names = side_chain_atoms
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        let roles = molframe::SideChainRoles {
            nitrogen,
            alpha_carbon,
            side_chain_atoms: &atom_names,
        };
        Ok(molframe::side_chain_definition(&component.0, &roles)
            .map(|definition| definition.atoms.iter().map(ToString::to_string).collect()))
    }

    fn peoe_charges(
        &self,
        py: Python<'_>,
        component_id: &str,
        options: PyPeoeOptions,
    ) -> PyResult<Vec<f64>> {
        let component = self
            .get(component_id)?
            .ok_or_else(|| PyKeyError::new_err(component_id.to_owned()))?;
        py.detach(move || molframe::component_peoe_charges(&component.0, options.0))
            .map_err(value_error)
    }
}

impl From<molframe::ComponentKind> for PyComponentKind {
    fn from(value: molframe::ComponentKind) -> Self {
        match value {
            molframe::ComponentKind::AminoAcid => Self::AminoAcid,
            molframe::ComponentKind::Nucleotide => Self::Nucleotide,
            molframe::ComponentKind::Saccharide => Self::Saccharide,
            molframe::ComponentKind::Lipid => Self::Lipid,
            molframe::ComponentKind::NonPolymer => Self::NonPolymer,
            molframe::ComponentKind::Solvent => Self::Solvent,
            molframe::ComponentKind::Ion => Self::Ion,
            molframe::ComponentKind::Unknown => Self::Unknown,
        }
    }
}

impl From<PyComponentKind> for molframe::ComponentKind {
    fn from(value: PyComponentKind) -> Self {
        match value {
            PyComponentKind::AminoAcid => Self::AminoAcid,
            PyComponentKind::Nucleotide => Self::Nucleotide,
            PyComponentKind::Saccharide => Self::Saccharide,
            PyComponentKind::Lipid => Self::Lipid,
            PyComponentKind::NonPolymer => Self::NonPolymer,
            PyComponentKind::Solvent => Self::Solvent,
            PyComponentKind::Ion => Self::Ion,
            PyComponentKind::Unknown => Self::Unknown,
        }
    }
}

impl From<molframe::StereoConfiguration> for PyStereoConfiguration {
    fn from(value: molframe::StereoConfiguration) -> Self {
        match value {
            molframe::StereoConfiguration::R => Self::R,
            molframe::StereoConfiguration::S => Self::S,
            molframe::StereoConfiguration::Mixed => Self::Mixed,
        }
    }
}

impl From<molframe::ComponentAtom> for PyComponentAtom {
    fn from(value: molframe::ComponentAtom) -> Self {
        Self {
            name: value.name.to_string(),
            alternate_name: value.alternate_name.map(|name| name.to_string()),
            element: PyElement(value.element),
            charge: value.charge,
            aromatic: value.aromatic,
            leaving: value.leaving,
            stereo: value.stereo.map(Into::into),
        }
    }
}

impl From<molframe::ComponentBond> for PyComponentBond {
    fn from(value: molframe::ComponentBond) -> Self {
        Self {
            atom_a: value.atom_a.to_string(),
            atom_b: value.atom_b.to_string(),
            order: value.order.into(),
            aromatic: value.aromatic,
            stereo: value.stereo.map(Into::into),
        }
    }
}

impl From<molframe::ComponentCoverage> for PyComponentCoverage {
    fn from(value: molframe::ComponentCoverage) -> Self {
        Self {
            coverage: value.coverage.into(),
            findings: value.findings.into_iter().map(Into::into).collect(),
            dictionary_version: value.dictionary_version.as_str().to_owned(),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyComponentKind>()?;
    module.add_class::<PyStereoConfiguration>()?;
    module.add_class::<PyComponentAtom>()?;
    module.add_class::<PyComponentBond>()?;
    module.add_class::<PyComponent>()?;
    module.add_class::<PyComponentCoverage>()?;
    module.add_class::<PyEquivalenceClasses>()?;
    module.add_class::<PyEquivalenceCache>()?;
    module.add_class::<PyPolymerAtomRole>()?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
