//! Actionable quality and backbone-validation records.

use crate::analysis::PyRamachandranOptions;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "QualityIssue", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyQualityIssue {
    ZeroOccupancy,
    OccupancyOutOfRange,
    NegativeBFactor,
}

#[pyclass(name = "QualityFlag", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyQualityFlag {
    atom: u32,
    issue: PyQualityIssue,
}

#[pymethods]
impl PyQualityFlag {
    #[getter]
    const fn atom(&self) -> u32 {
        self.atom
    }
    #[getter]
    const fn issue(&self) -> PyQualityIssue {
        self.issue
    }
}

#[pyclass(name = "RamachandranRegion", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyRamachandranRegion {
    AlphaHelixRight,
    BetaSheet,
    AlphaHelixLeft,
    Outlier,
}

#[pyclass(name = "RamachandranRecord", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRamachandranRecord {
    residue: u32,
    phi: f64,
    psi: f64,
    region: PyRamachandranRegion,
    reference_id: String,
    reference_version: String,
    distribution: String,
    probability: f64,
    percentile: Option<f64>,
}

#[pymethods]
impl PyRamachandranRecord {
    #[getter]
    const fn residue(&self) -> u32 {
        self.residue
    }
    #[getter]
    const fn phi(&self) -> f64 {
        self.phi
    }
    #[getter]
    const fn psi(&self) -> f64 {
        self.psi
    }
    #[getter]
    const fn region(&self) -> PyRamachandranRegion {
        self.region
    }
    #[getter]
    fn reference_id(&self) -> &str {
        &self.reference_id
    }
    #[getter]
    fn reference_version(&self) -> &str {
        &self.reference_version
    }
    #[getter]
    fn distribution(&self) -> &str {
        &self.distribution
    }
    #[getter]
    const fn probability(&self) -> f64 {
        self.probability
    }
    #[getter]
    const fn percentile(&self) -> Option<f64> {
        self.percentile
    }
}

impl From<molframe::validate::QualityFlag> for PyQualityFlag {
    fn from(value: molframe::validate::QualityFlag) -> Self {
        let issue = match value.issue {
            molframe::validate::QualityIssue::ZeroOccupancy => PyQualityIssue::ZeroOccupancy,
            molframe::validate::QualityIssue::OccupancyOutOfRange => {
                PyQualityIssue::OccupancyOutOfRange
            }
            molframe::validate::QualityIssue::NegativeBFactor => PyQualityIssue::NegativeBFactor,
        };
        Self {
            atom: value.atom.get(),
            issue,
        }
    }
}

impl PyRamachandranRecord {
    pub(super) fn rust_record(&self) -> molframe::validate::RamachandranRecord {
        molframe::validate::RamachandranRecord {
            residue: molframe::ResidueIndex::new(self.residue),
            phi: self.phi,
            psi: self.psi,
            region: match self.region {
                PyRamachandranRegion::AlphaHelixRight => {
                    molframe::validate::RamachandranRegion::AlphaHelixRight
                }
                PyRamachandranRegion::BetaSheet => {
                    molframe::validate::RamachandranRegion::BetaSheet
                }
                PyRamachandranRegion::AlphaHelixLeft => {
                    molframe::validate::RamachandranRegion::AlphaHelixLeft
                }
                PyRamachandranRegion::Outlier => molframe::validate::RamachandranRegion::Outlier,
            },
            reference: molframe::validate::ReferenceAssessment {
                set: self.reference_id.clone().into(),
                version: self.reference_version.clone().into(),
                distribution: self.distribution.clone().into(),
                probability: self.probability,
                percentile: self.percentile,
            },
        }
    }
}

impl From<molframe::validate::RamachandranRegion> for PyRamachandranRegion {
    fn from(value: molframe::validate::RamachandranRegion) -> Self {
        match value {
            molframe::validate::RamachandranRegion::AlphaHelixRight => Self::AlphaHelixRight,
            molframe::validate::RamachandranRegion::BetaSheet => Self::BetaSheet,
            molframe::validate::RamachandranRegion::AlphaHelixLeft => Self::AlphaHelixLeft,
            molframe::validate::RamachandranRegion::Outlier => Self::Outlier,
        }
    }
}

impl From<molframe::validate::RamachandranRecord> for PyRamachandranRecord {
    fn from(value: molframe::validate::RamachandranRecord) -> Self {
        Self {
            residue: value.residue.get(),
            phi: value.phi,
            psi: value.psi,
            region: value.region.into(),
            reference_id: value.reference.set.into(),
            reference_version: value.reference.version.into(),
            distribution: value.reference.distribution.into(),
            probability: value.reference.probability,
            percentile: value.reference.percentile,
        }
    }
}

#[pymethods]
impl PyStructure {
    fn quality_flags(&self, py: Python<'_>) -> Vec<PyQualityFlag> {
        let structure = self.structure().clone();
        py.detach(move || molframe::validate::quality_flags(&structure))
            .into_iter()
            .map(PyQualityFlag::from)
            .collect()
    }

    #[pyo3(signature = (options, *, outliers_only=false))]
    fn ramachandran(
        &self,
        py: Python<'_>,
        options: &PyRamachandranOptions,
        outliers_only: bool,
    ) -> PyResult<Vec<PyRamachandranRecord>> {
        let structure = self.structure().clone();
        let references = options.references.clone();
        let basins = options.basins.clone();
        let minimum_probability = options.minimum_probability;
        let values = py
            .detach(move || {
                let options = molframe::validate::RamachandranOptions::new(
                    &references,
                    basins,
                    minimum_probability,
                )
                .map_err(|error| error.to_string())?;
                if outliers_only {
                    molframe::validate::ramachandran_outliers(&structure, &options)
                } else {
                    molframe::validate::ramachandran(&structure, &options)
                }
                .map_err(|error| error.to_string())
            })
            .map_err(PyValueError::new_err)?;
        Ok(values
            .into_iter()
            .map(|value| PyRamachandranRecord {
                residue: value.residue.get(),
                phi: value.phi,
                psi: value.psi,
                region: match value.region {
                    molframe::validate::RamachandranRegion::AlphaHelixRight => {
                        PyRamachandranRegion::AlphaHelixRight
                    }
                    molframe::validate::RamachandranRegion::BetaSheet => {
                        PyRamachandranRegion::BetaSheet
                    }
                    molframe::validate::RamachandranRegion::AlphaHelixLeft => {
                        PyRamachandranRegion::AlphaHelixLeft
                    }
                    molframe::validate::RamachandranRegion::Outlier => {
                        PyRamachandranRegion::Outlier
                    }
                },
                reference_id: value.reference.set.into(),
                reference_version: value.reference.version.into(),
                distribution: value.reference.distribution.into(),
                probability: value.reference.probability,
                percentile: value.reference.percentile,
            })
            .collect())
    }
}
