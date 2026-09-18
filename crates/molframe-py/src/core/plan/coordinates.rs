//! Coordinate-backed declarative operations.

#[path = "coordinates/validation.rs"]
mod validation;

use super::dispatch::require_operation_tag;
use super::{Operation, PyPlan, execute_native};
use crate::contract::PyAnalysis;
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use validation::{
    contacts_complexity, validate_contacts, validate_coordinate_pair, validate_inclusion_radius,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum OperationKind {
    Contacts,
    Rmsd,
    Lddt,
    TmScore,
    GdtTs,
    GdtHa,
}

impl OperationKind {
    pub(super) const fn tag(self) -> &'static str {
        match self {
            Self::Contacts => "contacts",
            Self::Rmsd => "rmsd",
            Self::Lddt => "lddt",
            Self::TmScore => "tm_score",
            Self::GdtTs => "gdt_ts",
            Self::GdtHa => "gdt_ha",
        }
    }

    pub(super) fn parse(value: &str) -> Option<Self> {
        match value {
            "contacts" => Some(Self::Contacts),
            "rmsd" => Some(Self::Rmsd),
            "lddt" => Some(Self::Lddt),
            "tm_score" => Some(Self::TmScore),
            "gdt_ts" => Some(Self::GdtTs),
            "gdt_ha" => Some(Self::GdtHa),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum CoordinateMetric {
    Lddt { inclusion_radius: f64 },
    TmScore,
    GdtTs,
    GdtHa,
}

impl CoordinateMetric {
    fn operation(self) -> OperationKind {
        match self {
            Self::Lddt { .. } => OperationKind::Lddt,
            Self::TmScore => OperationKind::TmScore,
            Self::GdtTs => OperationKind::GdtTs,
            Self::GdtHa => OperationKind::GdtHa,
        }
    }

    pub(super) const fn native(self) -> molframe::ComparisonMetric {
        match self {
            Self::Lddt { inclusion_radius } => {
                molframe::ComparisonMetric::Lddt { inclusion_radius }
            }
            Self::TmScore => molframe::ComparisonMetric::TmScore,
            Self::GdtTs => molframe::ComparisonMetric::GdtTs,
            Self::GdtHa => molframe::ComparisonMetric::GdtHa,
        }
    }
}

pub(super) struct PyComparison {
    pub(super) mobile: Py<PyAny>,
    pub(super) reference: Py<PyAny>,
    pub(super) metric: CoordinateMetric,
}

#[pyclass(name = "Contacts", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyContacts {
    pub(crate) left: String,
    pub(crate) right: String,
    pub(crate) cutoff: f32,
    pub(crate) backend: PySpatialBackend,
    pub(crate) policy: PyAnalysisPolicy,
}

#[pymethods]
impl PyContacts {
    #[new]
    #[pyo3(signature = (*, left, right, cutoff=4.5, backend=None, policy=None))]
    fn new(
        left: String,
        right: String,
        cutoff: f32,
        backend: Option<PySpatialBackend>,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        validate_contacts(&left, &right, cutoff)?;
        Ok(Self {
            left,
            right,
            cutoff,
            backend: match backend {
                Some(backend) => backend,
                None => PySpatialBackend::Auto,
            },
            policy: match policy {
                Some(policy) => policy,
                None => PyAnalysisPolicy::new_default(),
            },
        })
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyAnalysis> {
        let plan = PyPlan {
            operations: vec![("result".to_owned(), Operation::Contacts(self.clone()))],
        };
        let result = execute_native(&plan, py, Some(structure), None)?;
        let Some(entry) = result.entries.into_iter().next() else {
            return Err(PyValueError::new_err(
                "native contacts plan returned no result",
            ));
        };
        match entry.value {
            molframe::PlanValue::Contacts(analysis) => {
                super::results::contact_analysis_to_py(py, *analysis)
            }
            _ => Err(PyValueError::new_err(
                "native contacts plan returned an incompatible result",
            )),
        }
    }

    pub(super) fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let result = PyDict::new(py);
        result.set_item("operation", OperationKind::Contacts.tag())?;
        result.set_item("backend", format!("{:?}", self.backend))?;
        result.set_item("complexity", contacts_complexity(self.backend))?;
        result.set_item("cutoff", self.cutoff)?;
        result.set_item("left", &self.left)?;
        result.set_item("right", &self.right)?;
        result.set_item("materializes", "ContactTable columns")?;
        result.set_item(
            "policy_fingerprint",
            self.policy.inner.fingerprint().to_string(),
        )?;
        Ok(result.unbind().into_any())
    }

    pub(super) fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let result = PyDict::new(py);
        result.set_item("operation", OperationKind::Contacts.tag())?;
        result.set_item("left", &self.left)?;
        result.set_item("right", &self.right)?;
        result.set_item("cutoff", self.cutoff)?;
        result.set_item("backend", self.backend)?;
        result.set_item("policy", Py::new(py, self.policy.clone())?)?;
        result.set_item(
            "policy_fingerprint",
            self.policy.inner.fingerprint().to_string(),
        )?;
        Ok(result.unbind().into_any())
    }

    #[staticmethod]
    pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, OperationKind::Contacts.tag())?;
        let left: String = super::required(config, "left")?.extract()?;
        let right: String = super::required(config, "right")?.extract()?;
        let cutoff: f32 = match super::optional(config, "cutoff")? {
            Some(value) => value.extract()?,
            None => 4.5,
        };
        let backend: PySpatialBackend = match super::optional(config, "backend")? {
            Some(value) => value.extract()?,
            None => PySpatialBackend::Auto,
        };
        let policy: PyAnalysisPolicy = match super::optional(config, "policy")? {
            Some(value) => value.extract()?,
            None => PyAnalysisPolicy::new_default(),
        };
        validate_contacts(&left, &right, cutoff)?;
        Ok(Self {
            left,
            right,
            cutoff,
            backend,
            policy,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Contacts(left={:?}, right={:?}, cutoff={}, backend={:?})",
            self.left, self.right, self.cutoff, self.backend
        )
    }
}

#[pyclass(name = "Rmsd", frozen)]
pub(crate) struct PyRmsd {
    pub(super) mobile: Py<PyAny>,
    pub(super) reference: Py<PyAny>,
}

#[pymethods]
impl PyRmsd {
    #[new]
    #[pyo3(signature = (*, mobile, reference))]
    fn new(py: Python<'_>, mobile: Py<PyAny>, reference: Py<PyAny>) -> PyResult<Self> {
        validate_coordinate_pair(py, &mobile, &reference)?;
        Ok(Self { mobile, reference })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<f64> {
        let plan = PyPlan {
            operations: vec![(
                "result".to_owned(),
                Operation::Rmsd(Self {
                    mobile: self.mobile.clone_ref(py),
                    reference: self.reference.clone_ref(py),
                }),
            )],
        };
        let result = execute_native(&plan, py, None, None)?;
        let Some(entry) = result.entries.into_iter().next() else {
            return Err(PyValueError::new_err("native RMSD plan returned no result"));
        };
        match entry.value {
            molframe::PlanValue::Rmsd(value) => Ok(value),
            _ => Err(PyValueError::new_err(
                "native RMSD plan returned an incompatible result",
            )),
        }
    }

    pub(super) fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let result = PyDict::new(py);
        result.set_item("operation", OperationKind::Rmsd.tag())?;
        result.set_item("backend", "rust")?;
        result.set_item("complexity", "O(n)")?;
        result.set_item("materializes", false)?;
        result.set_item("requires_c_contiguous", true)?;
        result.set_item("retained_input_count", 2)?;
        Ok(result.unbind().into_any())
    }

    pub(super) fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let result = PyDict::new(py);
        result.set_item("operation", OperationKind::Rmsd.tag())?;
        result.set_item("mobile", self.mobile.bind(py))?;
        result.set_item("reference", self.reference.bind(py))?;
        Ok(result.unbind().into_any())
    }

    #[staticmethod]
    pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, OperationKind::Rmsd.tag())?;
        Self::new(
            config.py(),
            super::required(config, "mobile")?.unbind(),
            super::required(config, "reference")?.unbind(),
        )
    }

    fn __repr__(&self) -> String {
        "Rmsd(inputs=2, mobile=<array>, reference=<array>)".to_owned()
    }
}

#[pyclass(name = "Lddt", frozen)]
pub(crate) struct PyLddt {
    pub(super) mobile: Py<PyAny>,
    pub(super) reference: Py<PyAny>,
    pub(super) inclusion_radius: f64,
}

#[pymethods]
impl PyLddt {
    #[new]
    #[pyo3(signature = (*, mobile, reference, inclusion_radius=15.0))]
    fn new(
        py: Python<'_>,
        mobile: Py<PyAny>,
        reference: Py<PyAny>,
        inclusion_radius: f64,
    ) -> PyResult<Self> {
        validate_inclusion_radius(inclusion_radius)?;
        validate_coordinate_pair(py, &mobile, &reference)?;
        Ok(Self {
            mobile,
            reference,
            inclusion_radius,
        })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<f64> {
        super::execute_comparison(
            py,
            &self.mobile,
            &self.reference,
            CoordinateMetric::Lddt {
                inclusion_radius: self.inclusion_radius,
            },
        )
    }

    pub(super) fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        comparison_explain(
            py,
            CoordinateMetric::Lddt {
                inclusion_radius: self.inclusion_radius,
            },
        )
    }

    pub(super) fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        comparison_to_dict(
            py,
            &self.mobile,
            &self.reference,
            CoordinateMetric::Lddt {
                inclusion_radius: self.inclusion_radius,
            },
        )
    }

    #[staticmethod]
    pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, OperationKind::Lddt.tag())?;
        let inclusion_radius: f64 = match super::optional(config, "inclusion_radius")? {
            Some(value) => value.extract()?,
            None => 15.0,
        };
        Self::new(
            config.py(),
            super::required(config, "mobile")?.unbind(),
            super::required(config, "reference")?.unbind(),
            inclusion_radius,
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "Lddt(inclusion_radius={}, mobile=<array>, reference=<array>)",
            self.inclusion_radius
        )
    }
}

macro_rules! simple_comparison_class {
    ($name:ident, $python:literal, $metric:expr) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $name {
            pub(super) mobile: Py<PyAny>,
            pub(super) reference: Py<PyAny>,
        }

        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (*, mobile, reference))]
            fn new(py: Python<'_>, mobile: Py<PyAny>, reference: Py<PyAny>) -> PyResult<Self> {
                validate_coordinate_pair(py, &mobile, &reference)?;
                Ok(Self { mobile, reference })
            }

            fn execute(&self, py: Python<'_>) -> PyResult<f64> {
                super::execute_comparison(py, &self.mobile, &self.reference, $metric)
            }

            pub(super) fn explain(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                comparison_explain(py, $metric)
            }

            pub(super) fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                comparison_to_dict(py, &self.mobile, &self.reference, $metric)
            }

            #[staticmethod]
            pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                require_operation_tag(config, ($metric).operation().tag())?;
                Self::new(
                    config.py(),
                    super::required(config, "mobile")?.unbind(),
                    super::required(config, "reference")?.unbind(),
                )
            }

            fn __repr__(&self) -> String {
                format!("{}(mobile=<array>, reference=<array>)", stringify!($name))
            }
        }
    };
}

simple_comparison_class!(PyTmScore, "TmScore", CoordinateMetric::TmScore);
simple_comparison_class!(PyGdtTs, "GdtTs", CoordinateMetric::GdtTs);
simple_comparison_class!(PyGdtHa, "GdtHa", CoordinateMetric::GdtHa);

pub(super) fn comparison_explain(py: Python<'_>, metric: CoordinateMetric) -> PyResult<Py<PyAny>> {
    let result = PyDict::new(py);
    result.set_item("operation", metric.operation().tag())?;
    result.set_item("backend", "rust")?;
    result.set_item(
        "complexity",
        if matches!(metric, CoordinateMetric::Lddt { .. }) {
            "O(n²)"
        } else {
            "O(n)"
        },
    )?;
    result.set_item("materializes", false)?;
    result.set_item("requires_c_contiguous", true)?;
    if let CoordinateMetric::Lddt { inclusion_radius } = metric {
        result.set_item("inclusion_radius", inclusion_radius)?;
    }
    Ok(result.unbind().into_any())
}

pub(super) fn comparison_to_dict(
    py: Python<'_>,
    mobile: &Py<PyAny>,
    reference: &Py<PyAny>,
    metric: CoordinateMetric,
) -> PyResult<Py<PyAny>> {
    let result = PyDict::new(py);
    result.set_item("operation", metric.operation().tag())?;
    result.set_item("mobile", mobile.bind(py))?;
    result.set_item("reference", reference.bind(py))?;
    if let CoordinateMetric::Lddt { inclusion_radius } = metric {
        result.set_item("inclusion_radius", inclusion_radius)?;
    }
    Ok(result.unbind().into_any())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyContacts>()?;
    module.add_class::<PyRmsd>()?;
    module.add_class::<PyLddt>()?;
    module.add_class::<PyTmScore>()?;
    module.add_class::<PyGdtTs>()?;
    module.add_class::<PyGdtHa>()?;
    Ok(())
}

#[cfg(test)]
#[path = "coordinates_tests.rs"]
mod tests;
