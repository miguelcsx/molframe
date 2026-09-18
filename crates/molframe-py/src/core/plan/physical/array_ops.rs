//! Declarative physical operations whose inputs are atom-aligned `NumPy` arrays.

use super::{PyPhysicalOperation, execute_operation};
use crate::analysis::{
    PyDensityGridSpec, PyLinearDensityOptions, PyPoreOptions, PySurfaceContactOptions,
};
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use numpy::PyReadonlyArray1;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

#[pyclass(name = "LinearDensity", frozen)]
pub(crate) struct PyLinearDensity {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyLinearDensity {
    #[new]
    #[pyo3(signature = (*, weights, options, policy=None))]
    fn new(
        py: Python<'_>,
        weights: Py<PyAny>,
        options: PyLinearDensityOptions,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        validate_scalar(py, &weights)?;
        Ok(Self {
            operation: PyPhysicalOperation::LinearDensity {
                weights,
                options: options.native(),
                policy: super::config::policy_or_default(policy),
            },
        })
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = super::config::from_dict(config)?;
        if let PyPhysicalOperation::LinearDensity { weights, .. } = &operation {
            validate_scalar(py, weights)?;
            Ok(Self { operation })
        } else {
            Err(super::config::expected("linear_density"))
        }
    }

    fn __repr__(&self) -> String {
        "LinearDensity(weights=<array>, options=<LinearDensityOptions>)".to_owned()
    }
}

#[pyclass(name = "_AnalysisDensityMap", frozen)]
pub(crate) struct PyAnalysisDensityMap {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyAnalysisDensityMap {
    #[new]
    #[pyo3(signature = (*, weights, spec, policy=None))]
    fn new(
        py: Python<'_>,
        weights: Py<PyAny>,
        spec: PyDensityGridSpec,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        validate_scalar(py, &weights)?;
        Ok(Self {
            operation: PyPhysicalOperation::DensityMap {
                weights,
                spec: spec.native(),
                policy: super::config::policy_or_default(policy),
            },
        })
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = super::config::from_dict(config)?;
        if let PyPhysicalOperation::DensityMap { weights, .. } = &operation {
            validate_scalar(py, weights)?;
            Ok(Self { operation })
        } else {
            Err(super::config::expected("density_map"))
        }
    }

    fn __repr__(&self) -> String {
        "DensityMap(weights=<array>, spec=<DensityGridSpec>)".to_owned()
    }
}

#[pyclass(name = "PoreProfile", frozen)]
pub(crate) struct PyPoreProfile {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PyPoreProfile {
    #[new]
    #[pyo3(signature = (*, radii, options, policy=None))]
    fn new(
        py: Python<'_>,
        radii: Py<PyAny>,
        options: PyPoreOptions,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        validate_float(py, &radii)?;
        Ok(Self {
            operation: PyPhysicalOperation::PoreProfile {
                radii,
                options: options.native(),
                policy: super::config::policy_or_default(policy),
            },
        })
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = super::config::from_dict(config)?;
        if let PyPhysicalOperation::PoreProfile { radii, .. } = &operation {
            validate_float(py, radii)?;
            Ok(Self { operation })
        } else {
            Err(super::config::expected("pore_profile"))
        }
    }

    fn __repr__(&self) -> String {
        "PoreProfile(radii=<array>, options=<PoreOptions>)".to_owned()
    }
}

#[pyclass(name = "SurfaceContacts", frozen)]
pub(crate) struct PySurfaceContacts {
    pub(crate) operation: PyPhysicalOperation,
}

#[pymethods]
impl PySurfaceContacts {
    #[new]
    #[pyo3(signature = (*, radii, options, policy=None))]
    fn new(
        py: Python<'_>,
        radii: Py<PyAny>,
        options: PySurfaceContactOptions,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        validate_float(py, &radii)?;
        let (tolerance, probe, density, minimum_area, backend) = options.native_parts();
        Ok(Self {
            operation: PyPhysicalOperation::SurfaceContacts {
                radii,
                tolerance,
                probe,
                density,
                minimum_area,
                backend,
                policy: super::config::policy_or_default(policy),
            },
        })
    }

    fn execute(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<crate::contract::PyAnalysis> {
        execute_operation(py, &self.operation, structure)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        super::config::operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = super::config::from_dict(config)?;
        if let PyPhysicalOperation::SurfaceContacts { radii, .. } = &operation {
            validate_float(py, radii)?;
            Ok(Self { operation })
        } else {
            Err(super::config::expected("surface_contacts"))
        }
    }

    fn __repr__(&self) -> String {
        "SurfaceContacts(radii=<array>, options=<SurfaceContactOptions>)".to_owned()
    }
}

fn validate_scalar(py: Python<'_>, value: &Py<PyAny>) -> PyResult<()> {
    let array = value.bind(py).extract::<PyReadonlyArray1<'_, f64>>()?;
    array.as_slice().map(|_| ()).map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "physical scalar inputs must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}

fn validate_float(py: Python<'_>, value: &Py<PyAny>) -> PyResult<()> {
    let array = value.bind(py).extract::<PyReadonlyArray1<'_, f32>>()?;
    array.as_slice().map(|_| ()).map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "physical float inputs must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}
