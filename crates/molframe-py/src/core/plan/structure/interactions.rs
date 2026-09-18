//! Declarative interaction operations.

use super::super::dispatch::require_operation_tag;
use super::{explain, from_dict, policy_or_default, run, to_dict};
use crate::analysis::{
    PyBasePairOptions, PyCationPiOptions, PyPiStackingOptions, PyWaterBridgeOptions,
};
use crate::chemistry::PyComponentDictionary;
use crate::contract::PyAnalysis;
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::sync::Arc;

#[pyclass(name = "BasePairs", frozen)]
pub(crate) struct PyBasePairs {
    pub(crate) request: molframe::StructureRequest,
    dictionary: Py<PyAny>,
}

#[pymethods]
impl PyBasePairs {
    #[new]
    #[pyo3(signature = (*, dictionary, options, policy=None))]
    fn new(
        py: Python<'_>,
        dictionary: Py<PyAny>,
        options: PyBasePairOptions,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        let provider = component_provider(py, &dictionary)?;
        Ok(Self {
            request: molframe::StructureRequest::BasePairs {
                provider,
                options: options.native(),
                policy: policy_or_default(policy),
            },
            dictionary,
        })
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
        run(py, &self.request, structure)
    }

    pub(crate) fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        explain(py, &self.request)
    }

    pub(crate) fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let config = to_dict(py, &self.request)?;
        let config = config.bind(py).cast::<PyDict>()?;
        config.set_item("dictionary", self.dictionary.bind(py))?;
        Ok(config.clone().unbind().into_any())
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let _ = py;
        from_dict_base_pairs(config)
    }

    fn __repr__(&self) -> String {
        "BasePairs()".to_owned()
    }
}

fn component_provider(
    py: Python<'_>,
    dictionary: &Py<PyAny>,
) -> PyResult<Arc<dyn molframe::ComponentProvider>> {
    component_provider_bound(dictionary.bind(py))
}

fn component_provider_bound(
    dictionary: &Bound<'_, PyAny>,
) -> PyResult<Arc<dyn molframe::ComponentProvider>> {
    let dictionary = dictionary.extract::<PyRef<'_, PyComponentDictionary>>()?;
    Ok(dictionary.0.clone())
}

#[pyclass(name = "HydrogenBonds", frozen)]
pub(crate) struct PyHydrogenBonds {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PyHydrogenBonds {
    #[new]
    #[pyo3(signature = (*, maximum_distance=3.5, minimum_angle=120.0, backend=None, periodic=false, policy=None))]
    fn new(
        maximum_distance: f32,
        minimum_angle: f64,
        backend: Option<PySpatialBackend>,
        periodic: bool,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::HydrogenBonds {
                options: molframe::analysis::HydrogenBondOptions {
                    maximum_donor_acceptor_distance: maximum_distance,
                    minimum_angle_degrees: minimum_angle,
                    backend: match backend {
                        Some(backend) => backend,
                        None => PySpatialBackend::Auto,
                    }
                    .into(),
                    periodic,
                },
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
            request: from_dict(config, "hydrogen_bonds")?,
        })
    }

    fn __repr__(&self) -> String {
        "HydrogenBonds()".to_owned()
    }
}

#[pyclass(name = "SaltBridges", frozen)]
pub(crate) struct PySaltBridges {
    pub(crate) request: molframe::StructureRequest,
}

#[pymethods]
impl PySaltBridges {
    #[new]
    #[pyo3(signature = (*, maximum_distance=4.0, backend=None, policy=None))]
    fn new(
        maximum_distance: f32,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            request: molframe::StructureRequest::SaltBridges {
                maximum_distance,
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
            request: from_dict(config, "salt_bridges")?,
        })
    }

    fn __repr__(&self) -> String {
        "SaltBridges()".to_owned()
    }
}

macro_rules! option_operation {
    ($name:ident, $python:literal, $tag:literal, $options:ty, $native:ident, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $name {
            pub(crate) request: molframe::StructureRequest,
        }

        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (*, options, policy=None))]
            fn new(options: $options, policy: Option<PyAnalysisPolicy>) -> Self {
                Self {
                    request: molframe::StructureRequest::$variant {
                        options: options.$native(),
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

option_operation!(
    PyPiStackingAnalysis,
    "PiStackingAnalysis",
    "pi_stacking",
    PyPiStackingOptions,
    native,
    PiStacking
);
option_operation!(
    PyCationPiAnalysis,
    "CationPiAnalysis",
    "cation_pi",
    PyCationPiOptions,
    native,
    CationPi
);
option_operation!(
    PyWaterBridges,
    "WaterBridges",
    "water_bridges",
    PyWaterBridgeOptions,
    native,
    WaterBridges
);

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
    extract!(PyHydrogenBonds);
    extract!(PySaltBridges);
    extract!(PyPiStackingAnalysis);
    extract!(PyCationPiAnalysis);
    extract!(PyWaterBridges);
    Ok(None)
}

pub(crate) fn typed_base_pairs(value: &Bound<'_, PyAny>) -> Option<PyBasePairs> {
    let Ok(operation) = value.extract::<PyRef<'_, PyBasePairs>>() else {
        return None;
    };
    Some(PyBasePairs {
        request: operation.request.clone(),
        dictionary: operation.dictionary.clone_ref(value.py()),
    })
}

pub(crate) fn serialized_base_pairs(
    py: Python<'_>,
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PyBasePairs>> {
    let Some(operation) = config.get_item("operation")? else {
        return Ok(None);
    };
    let operation: String = operation.extract()?;
    if operation != "base_pairs" {
        return Ok(None);
    }
    let _ = py;
    Ok(Some(from_dict_base_pairs(config)?))
}

fn from_dict_base_pairs(config: &Bound<'_, PyDict>) -> PyResult<PyBasePairs> {
    require_operation_tag(config, "base_pairs")?;
    let dictionary = config.get_item("dictionary")?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err("missing operation field: dictionary")
    })?;
    let provider = component_provider_bound(&dictionary)?;
    let request = super::config::from_dict_with_provider(config, Some(provider))?;
    Ok(PyBasePairs {
        request,
        dictionary: dictionary.unbind(),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBasePairs>()?;
    module.add_class::<PyHydrogenBonds>()?;
    module.add_class::<PySaltBridges>()?;
    module.add_class::<PyPiStackingAnalysis>()?;
    module.add_class::<PyCationPiAnalysis>()?;
    module.add_class::<PyWaterBridges>()?;
    Ok(())
}
