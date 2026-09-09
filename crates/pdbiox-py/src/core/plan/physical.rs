//! Declarative Python views over governed physical analyses.

use super::{Operation, PyPlan, execute_native};
use crate::analysis::{
    PyLeafletOptions, PyRadialOptions, coordination_analysis, density_map_analysis,
    leaflets_analysis, linear_density_analysis, pore_analysis, radial_analysis,
    surface_contacts_analysis,
};
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PySelection};
use crate::structure::PyStructure;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

#[path = "physical/array_ops.rs"]
mod array_ops;
#[path = "physical/config.rs"]
mod config;
#[path = "physical/native.rs"]
mod native;
pub(crate) use array_ops::{
    PyAnalysisDensityMap, PyLinearDensity, PyPoreProfile, PySurfaceContacts,
};
pub(super) use config::{explain_operation, operation_to_dict};
pub(crate) use config::{serialized_request, typed_request};
pub(crate) use native::to_native;

/// Owned native configuration for one physical analysis operation.
pub(crate) enum PyPhysicalOperation {
    RadialDistribution {
        left: pdbiox::AtomSelection,
        right: pdbiox::AtomSelection,
        options: pdbiox::analysis::RadialDistributionOptions,
        periodic: bool,
        policy: pdbiox::AnalysisPolicy,
    },
    CoordinationNumbers {
        left: pdbiox::AtomSelection,
        right: pdbiox::AtomSelection,
        minimum_distance: f32,
        maximum_distance: f32,
        backend: pdbiox::SpatialBackend,
        periodic: bool,
        policy: pdbiox::AnalysisPolicy,
    },
    Leaflets {
        sites: pdbiox::AtomSelection,
        options: pdbiox::analysis::LeafletOptions,
        periodic: bool,
        policy: pdbiox::AnalysisPolicy,
    },
    LinearDensity {
        weights: Py<PyAny>,
        options: pdbiox::analysis::LinearDensityOptions,
        policy: pdbiox::AnalysisPolicy,
    },
    DensityMap {
        weights: Py<PyAny>,
        spec: pdbiox::analysis::DensityGridSpec,
        policy: pdbiox::AnalysisPolicy,
    },
    PoreProfile {
        radii: Py<PyAny>,
        options: pdbiox::analysis::PoreProfileOptions,
        policy: pdbiox::AnalysisPolicy,
    },
    SurfaceContacts {
        radii: Py<PyAny>,
        tolerance: f32,
        probe: f32,
        density: f32,
        minimum_area: f32,
        backend: pdbiox::SpatialBackend,
        policy: pdbiox::AnalysisPolicy,
    },
}

impl PyPhysicalOperation {
    fn tag(&self) -> &'static str {
        match self {
            Self::RadialDistribution { .. } => "radial_distribution",
            Self::CoordinationNumbers { .. } => "coordination_numbers",
            Self::Leaflets { .. } => "leaflets",
            Self::LinearDensity { .. } => "linear_density",
            Self::DensityMap { .. } => "density_map",
            Self::PoreProfile { .. } => "pore_profile",
            Self::SurfaceContacts { .. } => "surface_contacts",
        }
    }

    fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::RadialDistribution {
                left,
                right,
                options,
                periodic,
                policy,
            } => Self::RadialDistribution {
                left: left.clone(),
                right: right.clone(),
                options: *options,
                periodic: *periodic,
                policy: policy.clone(),
            },
            Self::CoordinationNumbers {
                left,
                right,
                minimum_distance,
                maximum_distance,
                backend,
                periodic,
                policy,
            } => Self::CoordinationNumbers {
                left: left.clone(),
                right: right.clone(),
                minimum_distance: *minimum_distance,
                maximum_distance: *maximum_distance,
                backend: *backend,
                periodic: *periodic,
                policy: policy.clone(),
            },
            Self::Leaflets {
                sites,
                options,
                periodic,
                policy,
            } => Self::Leaflets {
                sites: sites.clone(),
                options: *options,
                periodic: *periodic,
                policy: policy.clone(),
            },
            Self::LinearDensity {
                weights,
                options,
                policy,
            } => Self::LinearDensity {
                weights: weights.clone_ref(py),
                options: *options,
                policy: policy.clone(),
            },
            Self::DensityMap {
                weights,
                spec,
                policy,
            } => Self::DensityMap {
                weights: weights.clone_ref(py),
                spec: *spec,
                policy: policy.clone(),
            },
            Self::PoreProfile {
                radii,
                options,
                policy,
            } => Self::PoreProfile {
                radii: radii.clone_ref(py),
                options: *options,
                policy: policy.clone(),
            },
            Self::SurfaceContacts {
                radii,
                tolerance,
                probe,
                density,
                minimum_area,
                backend,
                policy,
            } => Self::SurfaceContacts {
                radii: radii.clone_ref(py),
                tolerance: *tolerance,
                probe: *probe,
                density: *density,
                minimum_area: *minimum_area,
                backend: *backend,
                policy: policy.clone(),
            },
        }
    }
}

#[pyclass(name = "RadialDistribution", frozen)]
pub(crate) struct PyRadialDistribution {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyRadialDistribution {
    #[new]
    #[pyo3(signature = (*, left, right, options, periodic=false, policy=None))]
    fn new(
        left: &PySelection,
        right: &PySelection,
        options: PyRadialOptions,
        periodic: bool,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            operation: PyPhysicalOperation::RadialDistribution {
                left: left.inner.clone(),
                right: right.inner.clone(),
                options: options.native(),
                periodic,
                policy: config::policy_or_default(policy),
            },
        }
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = config::from_dict(config)?;
        match operation {
            PyPhysicalOperation::RadialDistribution { .. } => Ok(Self { operation }),
            _ => Err(config::expected("radial_distribution")),
        }
    }

    fn __repr__(&self) -> String {
        "RadialDistribution(left=<Selection>, right=<Selection>, options=<RadialOptions>)"
            .to_owned()
    }
}

#[pyclass(name = "CoordinationNumbers", frozen)]
pub(crate) struct PyCoordinationNumbers {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyCoordinationNumbers {
    #[new]
    #[pyo3(signature = (*, left, right, shell, backend=None, periodic=false, policy=None))]
    fn new(
        left: &PySelection,
        right: &PySelection,
        shell: (f32, f32),
        backend: Option<PySpatialBackend>,
        periodic: bool,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            operation: PyPhysicalOperation::CoordinationNumbers {
                left: left.inner.clone(),
                right: right.inner.clone(),
                minimum_distance: shell.0,
                maximum_distance: shell.1,
                backend: match backend {
                    Some(backend) => backend,
                    None => PySpatialBackend::Auto,
                }
                .into(),
                periodic,
                policy: config::policy_or_default(policy),
            },
        }
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = config::from_dict(config)?;
        match operation {
            PyPhysicalOperation::CoordinationNumbers { .. } => Ok(Self { operation }),
            _ => Err(config::expected("coordination_numbers")),
        }
    }

    fn __repr__(&self) -> String {
        "CoordinationNumbers(left=<Selection>, right=<Selection>, shell=...)".to_owned()
    }
}

#[pyclass(name = "Leaflets", frozen)]
pub(crate) struct PyLeaflets {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyLeaflets {
    #[new]
    #[pyo3(signature = (*, sites, options, periodic=false, policy=None))]
    fn new(
        sites: &PySelection,
        options: PyLeafletOptions,
        periodic: bool,
        policy: Option<PyAnalysisPolicy>,
    ) -> Self {
        Self {
            operation: PyPhysicalOperation::Leaflets {
                sites: sites.inner.clone(),
                options: options.native(),
                periodic,
                policy: config::policy_or_default(policy),
            },
        }
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = config::from_dict(config)?;
        match operation {
            PyPhysicalOperation::Leaflets { .. } => Ok(Self { operation }),
            _ => Err(config::expected("leaflets")),
        }
    }

    fn __repr__(&self) -> String {
        "Leaflets(sites=<Selection>, options=<LeafletOptions>)".to_owned()
    }
}

pub(super) fn execute_operation(
    py: Python<'_>,
    operation: &PyPhysicalOperation,
    structure: &PyStructure,
) -> PyResult<crate::contract::PyAnalysis> {
    let plan = PyPlan {
        operations: vec![(
            "result".to_owned(),
            Operation::Physical(operation.clone_ref(py)),
        )],
    };
    let result = execute_native(&plan, py, Some(structure), None)?;
    let Some(entry) = result.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native physical plan returned no result",
        ));
    };
    match entry.value {
        pdbiox::PlanValue::Physical(value) => value_to_python(py, *value),
        _ => Err(PyTypeError::new_err(
            "native physical plan returned an incompatible result",
        )),
    }
}

fn value_to_python(
    py: Python<'_>,
    value: pdbiox::PhysicalValue,
) -> PyResult<crate::contract::PyAnalysis> {
    match value {
        pdbiox::PhysicalValue::RadialDistribution(value) => radial_analysis(py, value),
        pdbiox::PhysicalValue::CoordinationNumbers(value) => coordination_analysis(py, value),
        pdbiox::PhysicalValue::Leaflets(value) => leaflets_analysis(py, value),
        pdbiox::PhysicalValue::LinearDensity(value) => linear_density_analysis(py, value),
        pdbiox::PhysicalValue::DensityMap(value) => density_map_analysis(py, value),
        pdbiox::PhysicalValue::PoreProfile(value) => pore_analysis(py, value),
        pdbiox::PhysicalValue::SurfaceContacts(value) => surface_contacts_analysis(py, value),
    }
}

pub(super) fn value_to_python_public(
    py: Python<'_>,
    value: pdbiox::PhysicalValue,
) -> PyResult<crate::contract::PyAnalysis> {
    value_to_python(py, value)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRadialDistribution>()?;
    module.add_class::<PyCoordinationNumbers>()?;
    module.add_class::<PyLeaflets>()?;
    module.add_class::<PyLinearDensity>()?;
    module.add_class::<PyAnalysisDensityMap>()?;
    module.add_class::<PyPoreProfile>()?;
    module.add_class::<PySurfaceContacts>()?;
    Ok(())
}
