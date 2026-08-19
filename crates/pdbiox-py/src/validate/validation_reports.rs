//! Lossless projections for validation reports that are not structure methods.

use super::{PyBondDeviation, PyNucleicTorsions, PyPucker};
use crate::geometry::PyEigenOptions;
use pyo3::prelude::*;

#[pyclass(name = "LigandGeometryReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLigandGeometryReport {
    #[pyo3(get)]
    pub(crate) outliers: Vec<PyBondDeviation>,
    #[pyo3(get)]
    pub(crate) intended: usize,
    #[pyo3(get)]
    pub(crate) assessed: usize,
}

#[pyclass(name = "RealSpaceCorrelation", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRealSpaceCorrelation {
    #[pyo3(get)]
    pub(crate) coefficient: f64,
    #[pyo3(get)]
    pub(crate) sample_count: usize,
    #[pyo3(get)]
    pub(crate) observed_mean: f64,
    #[pyo3(get)]
    pub(crate) calculated_mean: f64,
}

#[pyclass(name = "ReferenceGeometryOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyReferenceGeometryOptions(pub(crate) pdbiox::validate::ReferenceGeometryOptions);

#[pymethods]
impl PyReferenceGeometryOptions {
    #[new]
    fn new(maximum_bond_deviation: f64, maximum_angle_deviation_degrees: f64) -> Self {
        Self(pdbiox::validate::ReferenceGeometryOptions {
            maximum_bond_deviation,
            maximum_angle_deviation_degrees,
        })
    }

    #[getter]
    fn maximum_bond_deviation(&self) -> f64 {
        self.0.maximum_bond_deviation
    }

    #[getter]
    fn maximum_angle_deviation_degrees(&self) -> f64 {
        self.0.maximum_angle_deviation_degrees
    }
}

#[pyclass(name = "ReferenceBondFlag", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyReferenceBondFlag {
    #[pyo3(get)]
    pub(crate) residue: u32,
    #[pyo3(get)]
    pub(crate) first: u32,
    #[pyo3(get)]
    pub(crate) second: u32,
    #[pyo3(get)]
    pub(crate) observed: f64,
    #[pyo3(get)]
    pub(crate) expected: f64,
    #[pyo3(get)]
    pub(crate) deviation: f64,
}

#[pyclass(name = "ReferenceAngleFlag", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyReferenceAngleFlag {
    #[pyo3(get)]
    pub(crate) residue: u32,
    #[pyo3(get)]
    pub(crate) first: u32,
    #[pyo3(get)]
    pub(crate) centre: u32,
    #[pyo3(get)]
    pub(crate) third: u32,
    #[pyo3(get)]
    pub(crate) observed_degrees: f64,
    #[pyo3(get)]
    pub(crate) expected_degrees: f64,
    #[pyo3(get)]
    pub(crate) deviation_degrees: f64,
}

#[pyclass(name = "ReferenceGeometryReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceGeometryReport {
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyReferenceBondFlag>,
    #[pyo3(get)]
    pub(crate) angles: Vec<PyReferenceAngleFlag>,
    #[pyo3(get)]
    pub(crate) findings: Vec<String>,
    #[pyo3(get)]
    pub(crate) intended: usize,
    #[pyo3(get)]
    pub(crate) assessed: usize,
    #[pyo3(get)]
    pub(crate) options: PyReferenceGeometryOptions,
}

#[pyclass(name = "NucleicGeometryPolicy", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNucleicGeometryPolicy(pub(crate) pdbiox::validate::NucleicGeometryPolicy);

#[pymethods]
impl PyNucleicGeometryPolicy {
    #[new]
    fn new(
        maximum_base_plane_deviation: f64,
        glycosidic_bond_range: [f64; 2],
        phosphodiester_bond_range: [f64; 2],
        plane_fit: &PyEigenOptions,
    ) -> Self {
        Self(pdbiox::validate::NucleicGeometryPolicy {
            maximum_base_plane_deviation,
            glycosidic_bond_range,
            phosphodiester_bond_range,
            plane_fit: plane_fit.inner,
        })
    }

    #[getter]
    fn maximum_base_plane_deviation(&self) -> f64 {
        self.0.maximum_base_plane_deviation
    }

    #[getter]
    fn glycosidic_bond_range(&self) -> [f64; 2] {
        self.0.glycosidic_bond_range
    }

    #[getter]
    fn phosphodiester_bond_range(&self) -> [f64; 2] {
        self.0.phosphodiester_bond_range
    }

    #[getter]
    fn plane_fit(&self) -> PyEigenOptions {
        PyEigenOptions {
            inner: self.0.plane_fit,
        }
    }
}

#[pyclass(name = "NucleicGeometryIssue", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNucleicGeometryIssue {
    #[pyo3(get)]
    pub(crate) kind: String,
    #[pyo3(get)]
    pub(crate) available: Option<usize>,
    #[pyo3(get)]
    pub(crate) expected: Option<usize>,
    #[pyo3(get)]
    pub(crate) deviation: Option<f64>,
    #[pyo3(get)]
    pub(crate) distance: Option<f64>,
}

#[pyclass(name = "NucleicGeometryRecord", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNucleicGeometryRecord {
    #[pyo3(get)]
    pub(crate) residue: u32,
    #[pyo3(get)]
    pub(crate) torsions: PyNucleicTorsions,
    #[pyo3(get)]
    pub(crate) pucker: Option<PyPucker>,
    #[pyo3(get)]
    pub(crate) base_plane_deviation: Option<f64>,
    #[pyo3(get)]
    pub(crate) glycosidic_bond_length: Option<f64>,
    #[pyo3(get)]
    pub(crate) incoming_phosphodiester_length: Option<f64>,
    #[pyo3(get)]
    pub(crate) issues: Vec<PyNucleicGeometryIssue>,
}

impl From<pdbiox::validate::LigandGeometryReport> for PyLigandGeometryReport {
    fn from(value: pdbiox::validate::LigandGeometryReport) -> Self {
        Self {
            outliers: value.outliers.into_iter().map(Into::into).collect(),
            intended: value.intended,
            assessed: value.assessed,
        }
    }
}

impl From<pdbiox::validate::RealSpaceCorrelation> for PyRealSpaceCorrelation {
    fn from(value: pdbiox::validate::RealSpaceCorrelation) -> Self {
        Self {
            coefficient: value.coefficient,
            sample_count: value.sample_count,
            observed_mean: value.observed_mean,
            calculated_mean: value.calculated_mean,
        }
    }
}

impl From<pdbiox::validate::ReferenceBondFlag> for PyReferenceBondFlag {
    fn from(value: pdbiox::validate::ReferenceBondFlag) -> Self {
        Self {
            residue: value.residue.get(),
            first: value.first.get(),
            second: value.second.get(),
            observed: value.observed,
            expected: value.expected,
            deviation: value.deviation,
        }
    }
}

impl From<pdbiox::validate::ReferenceAngleFlag> for PyReferenceAngleFlag {
    fn from(value: pdbiox::validate::ReferenceAngleFlag) -> Self {
        Self {
            residue: value.residue.get(),
            first: value.first.get(),
            centre: value.centre.get(),
            third: value.third.get(),
            observed_degrees: value.observed_degrees,
            expected_degrees: value.expected_degrees,
            deviation_degrees: value.deviation_degrees,
        }
    }
}

impl From<pdbiox::validate::ReferenceGeometryReport> for PyReferenceGeometryReport {
    fn from(value: pdbiox::validate::ReferenceGeometryReport) -> Self {
        Self {
            bonds: value.bonds.into_iter().map(Into::into).collect(),
            angles: value.angles.into_iter().map(Into::into).collect(),
            findings: value
                .findings
                .into_iter()
                .map(|item| item.to_string())
                .collect(),
            intended: value.intended,
            assessed: value.assessed,
            options: PyReferenceGeometryOptions(value.options),
        }
    }
}

impl From<pdbiox::validate::NucleicGeometryIssue> for PyNucleicGeometryIssue {
    fn from(value: pdbiox::validate::NucleicGeometryIssue) -> Self {
        match value {
            pdbiox::validate::NucleicGeometryIssue::MissingBaseAtoms {
                available,
                expected,
            } => Self {
                kind: "missing_base_atoms".to_owned(),
                available: Some(available),
                expected: Some(expected),
                deviation: None,
                distance: None,
            },
            pdbiox::validate::NucleicGeometryIssue::MissingSugarAtoms { available } => Self {
                kind: "missing_sugar_atoms".to_owned(),
                available: Some(available),
                expected: None,
                deviation: None,
                distance: None,
            },
            pdbiox::validate::NucleicGeometryIssue::NonPlanarBase { deviation } => Self {
                kind: "non_planar_base".to_owned(),
                available: None,
                expected: None,
                deviation: Some(deviation),
                distance: None,
            },
            pdbiox::validate::NucleicGeometryIssue::GlycosidicBondLength { distance } => Self {
                kind: "glycosidic_bond_length".to_owned(),
                available: None,
                expected: None,
                deviation: None,
                distance: Some(distance),
            },
            pdbiox::validate::NucleicGeometryIssue::PhosphodiesterBondLength { distance } => Self {
                kind: "phosphodiester_bond_length".to_owned(),
                available: None,
                expected: None,
                deviation: None,
                distance: Some(distance),
            },
        }
    }
}

impl From<pdbiox::validate::NucleicGeometryRecord> for PyNucleicGeometryRecord {
    fn from(value: pdbiox::validate::NucleicGeometryRecord) -> Self {
        Self {
            residue: value.residue.get(),
            torsions: value.torsions.into(),
            pucker: value.pucker.map(Into::into),
            base_plane_deviation: value.base_plane_deviation,
            glycosidic_bond_length: value.glycosidic_bond_length,
            incoming_phosphodiester_length: value.incoming_phosphodiester_length,
            issues: value.issues.into_iter().map(Into::into).collect(),
        }
    }
}
