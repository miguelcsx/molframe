//! Declarative validation operations.

use super::{explain, from_dict, policy_or_default, run, to_dict};
use crate::analysis::PyPlanarityOptions;
use crate::chemistry::PyRadiusSet;
use crate::contract::PyAnalysis;
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

macro_rules! policy_operation {
    ($name:ident, $python:literal, $tag:literal, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $name {
            pub(crate) request: molframe::StructureRequest,
        }

        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (*, policy=None))]
            fn new(policy: Option<PyAnalysisPolicy>) -> Self {
                Self {
                    request: molframe::StructureRequest::$variant {
                        policy: policy_or_default(policy),
                    },
                }
            }

            fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
                run(py, &self.request, structure)
            }

            fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                explain(py, &self.request)
            }

            fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                to_dict(py, &self.request)
            }

            #[staticmethod]
            fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                Ok(Self {
                    request: from_dict(config, $tag)?,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}()", $python)
            }
        }
    };
}

policy_operation!(PyQualityFlags, "QualityFlags", "quality", Quality);
policy_operation!(PyValence, "Valence", "valence", Valence);
policy_operation!(PyCompleteness, "Completeness", "completeness", Completeness);

#[pyclass(name = "StericClashes", frozen)]
pub(crate) struct PyStericClashes {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyStericClashes {
    #[new]
    #[pyo3(signature = (*, tolerance=0.4, radii=None, backend=None, policy=None))]
    fn new(
        tolerance: f32,
        radii: Option<PyRadiusSet>,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::Clashes {
                tolerance,
                radii: match radii {
                    Some(radii) => radii,
                    None => PyRadiusSet::Bondi,
                }
                .into(),
                backend: match backend {
                    Some(backend) => backend,
                    None => PySpatialBackend::Auto,
                }
                .into(),
                policy: policy_or_default(policy),
            },
        }
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
        run(py, &self.request, structure)
    }

    fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        explain(py, &self.request)
    }

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_dict(py, &self.request)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        Ok(Self {
            request: from_dict(config, "clashes")?,
        })
    }

    fn __repr__(&self) -> String {
        "StericClashes()".to_owned()
    }
}

macro_rules! numeric_operation {
    ($name:ident, $python:literal, $tag:literal, $field:ident, $type:ty, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $name {
            pub(crate) request: molframe::StructureRequest,
        }

        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (*, $field, policy=None))]
            fn new($field: $type, policy: Option<PyAnalysisPolicy>) -> Self {
                Self {
                    request: molframe::StructureRequest::$variant {
                        $field,
                        policy: policy_or_default(policy),
                    },
                }
            }

            fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
                run(py, &self.request, structure)
            }

            fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                explain(py, &self.request)
            }

            fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                to_dict(py, &self.request)
            }

            #[staticmethod]
            fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                Ok(Self {
                    request: from_dict(config, $tag)?,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}()", $python)
            }
        }
    };
}

numeric_operation!(
    PyBondLengthDeviations,
    "BondLengthDeviations",
    "bond_length_deviations",
    tolerance,
    f32,
    BondLengthDeviations
);
numeric_operation!(
    PyCisPeptides,
    "CisPeptides",
    "cis_peptides",
    threshold_degrees,
    f64,
    CisPeptides
);

#[pyclass(name = "PlanarityCheck", frozen)]
pub(crate) struct PyPlanarityCheck {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyPlanarityCheck {
    #[new]
    #[pyo3(signature = (*, options, policy=None))]
    fn new(options: PyPlanarityOptions, policy: Option<PyAnalysisPolicy>) -> Self {
        Self {
            request: molframe::StructureRequest::Planarity {
                options: options.0,
                policy: policy_or_default(policy),
            },
        }
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
        run(py, &self.request, structure)
    }

    fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        explain(py, &self.request)
    }

    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_dict(py, &self.request)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        Ok(Self {
            request: from_dict(config, "planarity")?,
        })
    }

    fn __repr__(&self) -> String {
        "PlanarityCheck()".to_owned()
    }
}

pub(crate) fn typed_request(
    value: &Bound<'_, PyAny>,
) -> PyResult<Option<molframe::StructureRequest>> {
    macro_rules! extract {
        ($type:ty) => {
            if let Ok(operation) = value.extract::<PyRef<'_, $type>>() {
                return Ok(Some(operation.request.clone()));
            }
        };
    }
    extract!(PyQualityFlags);
    extract!(PyValence);
    extract!(PyCompleteness);
    extract!(PyStericClashes);
    extract!(PyBondLengthDeviations);
    extract!(PyCisPeptides);
    extract!(PyPlanarityCheck);
    Ok(None)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyQualityFlags>()?;
    module.add_class::<PyValence>()?;
    module.add_class::<PyCompleteness>()?;
    module.add_class::<PyStericClashes>()?;
    module.add_class::<PyBondLengthDeviations>()?;
    module.add_class::<PyCisPeptides>()?;
    module.add_class::<PyPlanarityCheck>()?;
    Ok(())
}
