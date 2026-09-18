//! Mechanical Python projections of native geometric validation results.

use crate::chemistry::PyRadiusSet;
use crate::geometry::PyEigenOptions;
use crate::graph::PySpatialBackend;
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::{PyRamachandranRecord, PyReferenceLibrary};

#[pyclass(name = "BondDeviation", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBondDeviation {
    inner: molframe::validate::BondDeviation,
    #[pyo3(get)]
    atom_a: u32,
    #[pyo3(get)]
    atom_b: u32,
    #[pyo3(get)]
    observed: f32,
    #[pyo3(get)]
    expected: f32,
    #[pyo3(get)]
    deviation: f32,
}

#[pyclass(name = "Clash", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyClash {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
    #[pyo3(get)]
    overlap: f32,
}

#[pyclass(name = "MissingResidue", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMissingResidue {
    #[pyo3(get)]
    canonical_position: u32,
    #[pyo3(get)]
    component: String,
}

impl PyMissingResidue {
    pub(crate) fn from_sequence(value: molframe::MissingResidue, component: &str) -> Self {
        Self {
            canonical_position: value.canonical_position,
            component: component.to_owned(),
        }
    }
}

#[pyclass(name = "ChainCompleteness", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainCompleteness {
    #[pyo3(get)]
    chain: String,
    #[pyo3(get)]
    observed: usize,
    #[pyo3(get)]
    canonical: usize,
    #[pyo3(get)]
    missing: Vec<PyMissingResidue>,
}

#[pyclass(name = "CisPeptide", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCisPeptide {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    omega: f64,
}

#[pyclass(name = "PlanarityFlag", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPlanarityFlag {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    deviation: f64,
}

#[pyclass(name = "PlanarityOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPlanarityOptions(pub(crate) molframe::validate::PlanarityOptions);

#[pymethods]
impl PyPlanarityOptions {
    #[new]
    fn new(maximum_deviation: f64, plane_fit: &PyEigenOptions) -> Self {
        Self(molframe::validate::PlanarityOptions {
            maximum_deviation,
            plane_fit: plane_fit.inner,
        })
    }
}

impl PyPlanarityOptions {
    pub(crate) const fn from_native_parts(value: molframe::validate::PlanarityOptions) -> Self {
        Self(value)
    }
}

#[pyclass(name = "ValenceError", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyValenceError {
    #[pyo3(get)]
    atom: u32,
    #[pyo3(get)]
    bonds: usize,
    #[pyo3(get)]
    maximum: usize,
}

#[pyclass(name = "ReferenceAssessment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceAssessment {
    #[pyo3(get)]
    set: String,
    #[pyo3(get)]
    version: String,
    #[pyo3(get)]
    distribution: String,
    #[pyo3(get)]
    probability: f64,
    #[pyo3(get)]
    percentile: Option<f64>,
}

#[pymethods]
impl PyStructure {
    fn bond_length_deviations(&self, py: Python<'_>, tolerance: f32) -> Vec<PyBondDeviation> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::bond_length_deviations(&structure, tolerance))
            .into_iter()
            .map(PyBondDeviation::from)
            .collect()
    }

    fn ligand_geometry_outliers(&self, py: Python<'_>, tolerance: f32) -> Vec<PyBondDeviation> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::ligand_geometry_outliers(&structure, tolerance))
            .into_iter()
            .map(PyBondDeviation::from)
            .collect()
    }

    #[pyo3(signature = (tolerance, radius_set, *, backend=PySpatialBackend::Auto, context=None))]
    fn clashes(
        &self,
        py: Python<'_>,
        tolerance: f32,
        radius_set: PyRadiusSet,
        backend: PySpatialBackend,
        context: Option<&crate::core::execution::PyExecutionContext>,
    ) -> PyResult<Vec<PyClash>> {
        let structure = self.structure().clone();
        let context = context.map_or_else(
            crate::core::execution::default_context,
            crate::core::execution::PyExecutionContext::native,
        );
        py.detach(move || {
            molframe::validate::clashes(
                &structure,
                tolerance,
                radius_set.into(),
                backend.into(),
                &context,
            )
        })
        .map(|values| values.into_iter().map(PyClash::from).collect())
        .map_err(value_error)
    }

    fn completeness(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
    ) -> PyResult<Vec<PyChainCompleteness>> {
        let structure = self.structure().clone();
        py.detach(move || {
            molframe::validate::completeness(&structure, namespace.into())
                .map_err(|error| error.to_string())?
                .into_iter()
                .map(|value| {
                    let missing = value
                        .missing
                        .into_iter()
                        .map(|missing| {
                            structure
                                .resolve(missing.component)
                                .map(|component| PyMissingResidue {
                                    canonical_position: missing.canonical_position,
                                    component: component.to_owned(),
                                })
                                .ok_or_else(|| {
                                    "canonical component symbol is not present in the structure"
                                        .to_owned()
                                })
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    Ok(PyChainCompleteness {
                        chain: value.chain,
                        observed: value.observed,
                        canonical: value.canonical,
                        missing,
                    })
                })
                .collect::<Result<Vec<_>, String>>()
        })
        .map_err(value_error)
    }

    fn cis_peptides(&self, py: Python<'_>, threshold_degrees: f64) -> PyResult<Vec<PyCisPeptide>> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::cis_peptides(&structure, threshold_degrees))
            .map(|values| values.into_iter().map(PyCisPeptide::from).collect())
            .map_err(value_error)
    }

    fn nonplanar_aromatic_rings(
        &self,
        py: Python<'_>,
        options: PyPlanarityOptions,
    ) -> PyResult<Vec<PyPlanarityFlag>> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::nonplanar_aromatic_rings(&structure, options.0))
            .map(|values| values.into_iter().map(PyPlanarityFlag::from).collect())
            .map_err(value_error)
    }

    fn overvalent_atoms(&self, py: Python<'_>) -> Vec<PyValenceError> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::overvalent_atoms(&structure))
            .into_iter()
            .map(PyValenceError::from)
            .collect()
    }
}

#[pyfunction]
pub(crate) fn assess_bond_deviation(
    py: Python<'_>,
    deviation: &PyBondDeviation,
    references: &PyReferenceLibrary,
    distribution: &str,
) -> PyResult<PyReferenceAssessment> {
    let deviation = deviation.inner;
    let references = references.0.clone();
    let distribution = distribution.to_owned();
    py.detach(move || {
        molframe::validate::assess_bond_deviation(&deviation, &references, &distribution)
    })
    .map(PyReferenceAssessment::from)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn assess_ramachandran(
    py: Python<'_>,
    record: &PyRamachandranRecord,
    references: &PyReferenceLibrary,
    distribution: &str,
) -> PyResult<PyReferenceAssessment> {
    let record = record.rust_record();
    let references = references.0.clone();
    let distribution = distribution.to_owned();
    py.detach(move || molframe::validate::assess_ramachandran(&record, &references, &distribution))
        .map(PyReferenceAssessment::from)
        .map_err(value_error)
}

macro_rules! project {
    ($rust:ty, $python:ty, $body:expr) => {
        impl From<$rust> for $python {
            fn from(value: $rust) -> Self {
                $body(value)
            }
        }
    };
}

project!(
    molframe::validate::BondDeviation,
    PyBondDeviation,
    |value: molframe::validate::BondDeviation| {
        Self {
            inner: value,
            atom_a: value.atom_a.get(),
            atom_b: value.atom_b.get(),
            observed: value.observed,
            expected: value.expected,
            deviation: value.deviation,
        }
    }
);
project!(
    molframe::validate::Clash,
    PyClash,
    |value: molframe::validate::Clash| Self {
        first: value.first.get(),
        second: value.second.get(),
        overlap: value.overlap,
    }
);
project!(
    molframe::validate::CisPeptide,
    PyCisPeptide,
    |value: molframe::validate::CisPeptide| Self {
        residue: value.residue.get(),
        omega: value.omega,
    }
);
project!(
    molframe::validate::PlanarityFlag,
    PyPlanarityFlag,
    |value: molframe::validate::PlanarityFlag| {
        Self {
            residue: value.residue.get(),
            deviation: value.deviation,
        }
    }
);
project!(
    molframe::validate::ValenceError,
    PyValenceError,
    |value: molframe::validate::ValenceError| {
        Self {
            atom: value.atom.get(),
            bonds: value.bonds,
            maximum: value.maximum,
        }
    }
);
project!(
    molframe::validate::ReferenceAssessment,
    PyReferenceAssessment,
    |value: molframe::validate::ReferenceAssessment| Self {
        set: value.set.into(),
        version: value.version.into(),
        distribution: value.distribution.into(),
        probability: value.probability,
        percentile: value.percentile,
    }
);

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn project_completeness(
    structure: &molframe::Structure,
    value: molframe::validate::ChainCompleteness,
) -> Result<PyChainCompleteness, String> {
    let missing = value
        .missing
        .into_iter()
        .map(|missing| {
            structure
                .resolve(missing.component)
                .map(|component| PyMissingResidue {
                    canonical_position: missing.canonical_position,
                    component: component.to_owned(),
                })
                .ok_or_else(|| {
                    "canonical component symbol is not present in the structure".to_owned()
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(PyChainCompleteness {
        chain: value.chain,
        observed: value.observed,
        canonical: value.canonical,
        missing,
    })
}
