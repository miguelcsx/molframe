//! Owned Python projections for typed `ModelCIF` metadata and QA metrics.

use crate::cif_document::{PyCifDocument, PyCifValue};
use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "ModelRow", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyModelRow {
    #[pyo3(get)]
    values: Vec<PyCifValue>,
}

#[pymethods]
impl PyModelRow {
    fn value(&self, index: usize) -> Option<PyCifValue> {
        self.values.get(index).cloned()
    }
}

#[pyclass(name = "ModelCategory", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyModelCategory {
    inner: pdbiox::ModelCategory,
}

#[pymethods]
impl PyModelCategory {
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn items(&self) -> Vec<String> {
        self.inner.items().iter().map(ToString::to_string).collect()
    }

    #[getter]
    fn rows(&self) -> Vec<PyModelRow> {
        self.inner
            .rows()
            .iter()
            .map(|row| PyModelRow {
                values: self
                    .inner
                    .items()
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        row.value(index).cloned().map_or_else(
                            || PyCifValue {
                                inner: pdbiox::CifValue::Unknown,
                            },
                            |inner| PyCifValue { inner },
                        )
                    })
                    .collect(),
            })
            .collect()
    }

    fn value(&self, item: &str, row: usize) -> Option<PyCifValue> {
        self.inner
            .value(item, row)
            .cloned()
            .map(|inner| PyCifValue { inner })
    }
}

#[pyclass(name = "MetricDefinition", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMetricDefinition {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    name: Option<String>,
    #[pyo3(get)]
    metric_type: String,
    #[pyo3(get)]
    mode: String,
    #[pyo3(get)]
    software_group_id: Option<String>,
}

#[pyclass(name = "GlobalMetric", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGlobalMetric {
    #[pyo3(get)]
    model_id: String,
    #[pyo3(get)]
    metric_id: String,
    #[pyo3(get)]
    value: f64,
}

#[pyclass(name = "LocalMetric", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLocalMetric {
    #[pyo3(get)]
    model_id: String,
    #[pyo3(get)]
    chain_id: String,
    #[pyo3(get)]
    sequence_id: i64,
    #[pyo3(get)]
    component_id: Option<String>,
    #[pyo3(get)]
    metric_id: String,
    #[pyo3(get)]
    value: f64,
}

#[pyclass(name = "PairwiseMetric", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPairwiseMetric {
    #[pyo3(get)]
    model_id: String,
    #[pyo3(get)]
    first_chain_id: String,
    #[pyo3(get)]
    first_sequence_id: i64,
    #[pyo3(get)]
    second_chain_id: String,
    #[pyo3(get)]
    second_sequence_id: i64,
    #[pyo3(get)]
    metric_id: String,
    #[pyo3(get)]
    value: f64,
}

#[pyclass(name = "QualityMetrics", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyQualityMetrics {
    inner: pdbiox::QualityMetrics,
}

#[pymethods]
impl PyQualityMetrics {
    #[getter]
    fn definitions(&self) -> Vec<PyMetricDefinition> {
        self.inner
            .definitions
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn global_metrics(&self) -> Vec<PyGlobalMetric> {
        self.inner.global.iter().cloned().map(Into::into).collect()
    }

    #[getter]
    fn local(&self) -> Vec<PyLocalMetric> {
        self.inner.local.iter().cloned().map(Into::into).collect()
    }

    #[getter]
    fn pairwise(&self) -> Vec<PyPairwiseMetric> {
        self.inner
            .pairwise
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn plddt(&self) -> Vec<PyLocalMetric> {
        self.inner.plddt().cloned().map(Into::into).collect()
    }

    fn pae(&self) -> Vec<PyPairwiseMetric> {
        self.inner.pae().cloned().map(Into::into).collect()
    }

    fn ptm(&self) -> Vec<PyGlobalMetric> {
        self.inner.ptm().cloned().map(Into::into).collect()
    }
}

macro_rules! metadata_projection {
    ($name:ident, $python:literal, $native:ty { $( $field:ident : $type:ty ),+ $(,)? }) => {
        #[pyclass(name = $python, frozen, skip_from_py_object)]
        #[derive(Clone, Debug)]
        pub(crate) struct $name {
            $(#[pyo3(get)] $field: $type,)+
        }
    };
}

metadata_projection!(PyModelDescription, "ModelDescription", pdbiox::ModelDescription {
    id: String,
    assembly_id: Option<String>,
    name: Option<String>,
    model_type: Option<String>,
});
metadata_projection!(PyTarget, "Target", pdbiox::Target {
    entity_id: String,
    data_id: Option<String>,
    origin: Option<String>,
});
metadata_projection!(PyTemplate, "Template", pdbiox::Template {
    id: String,
    target_chain_id: Option<String>,
    name: Option<String>,
    origin: Option<String>,
    entity_type: Option<String>,
});
metadata_projection!(PyProtocolStep, "ProtocolStep", pdbiox::ProtocolStep {
    protocol_id: String,
    step_id: String,
    method_type: String,
    name: Option<String>,
    details: Option<String>,
    software_group_id: Option<String>,
});

#[pyclass(name = "SoftwareGroup", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySoftwareGroup {
    group: String,
    software: String,
    parameter_group: Option<String>,
}

#[pymethods]
impl PySoftwareGroup {
    #[getter]
    fn group_id(&self) -> &str {
        &self.group
    }

    #[getter]
    fn software_id(&self) -> &str {
        &self.software
    }

    #[getter]
    fn parameter_group_id(&self) -> Option<&str> {
        self.parameter_group.as_deref()
    }
}

#[pyclass(name = "ModelCif", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyModelCif {
    pub(crate) inner: pdbiox::ModelCif,
}

impl PyModelCif {
    pub(crate) fn inner(&self) -> &pdbiox::ModelCif {
        &self.inner
    }
}

#[pymethods]
impl PyModelCif {
    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[getter]
    fn categories(&self) -> Vec<PyModelCategory> {
        self.inner
            .categories
            .iter()
            .cloned()
            .map(|inner| PyModelCategory { inner })
            .collect()
    }

    #[getter]
    fn models(&self) -> Vec<PyModelDescription> {
        self.inner.models.iter().cloned().map(Into::into).collect()
    }

    #[getter]
    fn targets(&self) -> Vec<PyTarget> {
        self.inner.targets.iter().cloned().map(Into::into).collect()
    }

    #[getter]
    fn templates(&self) -> Vec<PyTemplate> {
        self.inner
            .templates
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn protocol_steps(&self) -> Vec<PyProtocolStep> {
        self.inner
            .protocol_steps
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn software_groups(&self) -> Vec<PySoftwareGroup> {
        self.inner
            .software_groups
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn quality(&self) -> PyQualityMetrics {
        PyQualityMetrics {
            inner: self.inner.quality.clone(),
        }
    }

    fn confidence(&self) -> PyQualityMetrics {
        PyQualityMetrics::from_inner(self.inner.quality.clone())
    }
}

#[pymethods]
impl PyStructure {
    fn model_cif(&self) -> Option<PyModelCif> {
        self.structure()
            .extensions()
            .get::<pdbiox::ModelCif>(pdbiox::MODEL_CIF_EXTENSION)
            .cloned()
            .map(|inner| PyModelCif { inner })
    }

    fn confidence(&self) -> Option<PyQualityMetrics> {
        self.model_cif()
            .map(|model| PyQualityMetrics::from_inner(model.inner.quality))
    }
}

#[pyfunction]
pub(crate) fn lower_model_cif(py: Python<'_>, document: &PyCifDocument) -> PyResult<PyModelCif> {
    let (model, findings) = pdbiox::modelcif::lower(&document.inner);
    if findings.is_empty() {
        Ok(PyModelCif { inner: model })
    } else {
        Err(read_error(py, &findings))
    }
}

impl PyQualityMetrics {
    fn from_inner(inner: pdbiox::QualityMetrics) -> Self {
        Self { inner }
    }
}

impl From<pdbiox::MetricDefinition> for PyMetricDefinition {
    fn from(value: pdbiox::MetricDefinition) -> Self {
        Self {
            id: value.id.to_string(),
            name: value.name.map(|value| value.to_string()),
            metric_type: value.metric_type.to_string(),
            mode: value.mode.to_string(),
            software_group_id: value.software_group_id.map(|value| value.to_string()),
        }
    }
}

impl From<pdbiox::GlobalMetric> for PyGlobalMetric {
    fn from(value: pdbiox::GlobalMetric) -> Self {
        Self {
            model_id: value.model_id.to_string(),
            metric_id: value.metric_id.to_string(),
            value: value.value,
        }
    }
}

impl From<pdbiox::LocalMetric> for PyLocalMetric {
    fn from(value: pdbiox::LocalMetric) -> Self {
        Self {
            model_id: value.model_id.to_string(),
            chain_id: value.chain_id.to_string(),
            sequence_id: value.sequence_id,
            component_id: value.component_id.map(|value| value.to_string()),
            metric_id: value.metric_id.to_string(),
            value: value.value,
        }
    }
}

impl From<pdbiox::PairwiseMetric> for PyPairwiseMetric {
    fn from(value: pdbiox::PairwiseMetric) -> Self {
        Self {
            model_id: value.model_id.to_string(),
            first_chain_id: value.first_chain_id.to_string(),
            first_sequence_id: value.first_sequence_id,
            second_chain_id: value.second_chain_id.to_string(),
            second_sequence_id: value.second_sequence_id,
            metric_id: value.metric_id.to_string(),
            value: value.value,
        }
    }
}

impl From<pdbiox::ModelDescription> for PyModelDescription {
    fn from(value: pdbiox::ModelDescription) -> Self {
        Self {
            id: value.id.to_string(),
            assembly_id: value.assembly_id.map(|value| value.to_string()),
            name: value.name.map(|value| value.to_string()),
            model_type: value.model_type.map(|value| value.to_string()),
        }
    }
}

impl From<pdbiox::Target> for PyTarget {
    fn from(value: pdbiox::Target) -> Self {
        Self {
            entity_id: value.entity_id.to_string(),
            data_id: value.data_id.map(|value| value.to_string()),
            origin: value.origin.map(|value| value.to_string()),
        }
    }
}

impl From<pdbiox::Template> for PyTemplate {
    fn from(value: pdbiox::Template) -> Self {
        Self {
            id: value.id.to_string(),
            target_chain_id: value.target_chain_id.map(|value| value.to_string()),
            name: value.name.map(|value| value.to_string()),
            origin: value.origin.map(|value| value.to_string()),
            entity_type: value.entity_type.map(|value| value.to_string()),
        }
    }
}

impl From<pdbiox::ProtocolStep> for PyProtocolStep {
    fn from(value: pdbiox::ProtocolStep) -> Self {
        Self {
            protocol_id: value.protocol_id.to_string(),
            step_id: value.step_id.to_string(),
            method_type: value.method_type.to_string(),
            name: value.name.map(|value| value.to_string()),
            details: value.details.map(|value| value.to_string()),
            software_group_id: value.software_group_id.map(|value| value.to_string()),
        }
    }
}

impl From<pdbiox::SoftwareGroup> for PySoftwareGroup {
    fn from(value: pdbiox::SoftwareGroup) -> Self {
        Self {
            group: value.group_id.to_string(),
            software: value.software_id.to_string(),
            parameter_group: value.parameter_group_id.map(|value| value.to_string()),
        }
    }
}
