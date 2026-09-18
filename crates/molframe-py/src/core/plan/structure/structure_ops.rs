//! Declarative structure-level analysis operations.

use super::{explain, from_dict, policy_or_default, run, to_dict};
use crate::analysis::{PyDsspOptions, PyGnmOptions};
use crate::contract::PyAnalysis;
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PySelection};
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

policy_operation!(
    PyNucleicTorsionAnalysis,
    "NucleicTorsionAnalysis",
    "nucleic_torsions",
    NucleicTorsions
);

#[pyclass(name = "ResidueContacts", frozen)]
pub(crate) struct PyResidueContacts {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyResidueContacts {
    #[new]
    #[pyo3(signature = (*, cutoff=4.5, minimum_separation=0, backend=None, policy=None))]
    fn new(
        cutoff: f32,
        minimum_separation: u32,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::ContactMap {
                cutoff,
                minimum_separation,
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
            request: from_dict(config, "contact_map")?,
        })
    }

    fn __repr__(&self) -> String {
        "ResidueContacts()".to_owned()
    }
}

#[pyclass(name = "ChainInterfaceAnalysis", frozen)]
pub(crate) struct PyChainInterfaceAnalysis {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyChainInterfaceAnalysis {
    #[new]
    #[pyo3(signature = (*, first_chain, second_chain, cutoff=4.5, backend=None, policy=None))]
    fn new(
        first_chain: String,
        second_chain: String,
        cutoff: f32,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::ChainInterface {
                first_chain: first_chain.into_boxed_str(),
                second_chain: second_chain.into_boxed_str(),
                cutoff,
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
            request: from_dict(config, "chain_interface")?,
        })
    }

    fn __repr__(&self) -> String {
        "ChainInterfaceAnalysis()".to_owned()
    }
}

#[pyclass(name = "AssignSecondaryStructure", frozen)]
pub(crate) struct PyAssignSecondaryStructure {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyAssignSecondaryStructure {
    #[new]
    #[pyo3(signature = (*, options, policy=None))]
    fn new(options: PyDsspOptions, policy: Option<PyAnalysisPolicy>) -> Self {
        Self {
            request: molframe::StructureRequest::SecondaryStructure {
                options: options.native(),
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
            request: from_dict(config, "secondary_structure")?,
        })
    }

    fn __repr__(&self) -> String {
        "AssignSecondaryStructure()".to_owned()
    }
}

#[pyclass(name = "HalfSphereExposureAnalysis", frozen)]
pub(crate) struct PyHalfSphereExposureAnalysis {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyHalfSphereExposureAnalysis {
    #[new]
    #[pyo3(signature = (*, radius, backend=None, policy=None))]
    fn new(
        radius: f32,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::HalfSphereExposure {
                radius,
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
            request: from_dict(config, "half_sphere_exposure")?,
        })
    }

    fn __repr__(&self) -> String {
        "HalfSphereExposureAnalysis()".to_owned()
    }
}

#[pyclass(name = "Gnm", frozen)]
pub(crate) struct PyGnm {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyGnm {
    #[new]
    #[pyo3(signature = (*, sites, options, periodic=false, policy=None))]
    fn new(
        sites: &PySelection,
        options: PyGnmOptions,
        periodic: bool,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::GaussianNetworkModel {
                sites: sites.inner.clone(),
                options: options.native(),
                periodic,
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
            request: from_dict(config, "gaussian_network_model")?,
        })
    }

    fn __repr__(&self) -> String {
        "Gnm()".to_owned()
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
    extract!(PyNucleicTorsionAnalysis);
    extract!(PyResidueContacts);
    extract!(PyChainInterfaceAnalysis);
    extract!(PyAssignSecondaryStructure);
    extract!(PyHalfSphereExposureAnalysis);
    extract!(PyGnm);
    Ok(None)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNucleicTorsionAnalysis>()?;
    module.add_class::<PyResidueContacts>()?;
    module.add_class::<PyChainInterfaceAnalysis>()?;
    module.add_class::<PyAssignSecondaryStructure>()?;
    module.add_class::<PyHalfSphereExposureAnalysis>()?;
    module.add_class::<PyGnm>()?;
    Ok(())
}
