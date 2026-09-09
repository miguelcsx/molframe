//! Owned Python projections for typed `ModelCIF` metadata and QA metrics.

use crate::cif_document::{PyCifDocument, PyCifValue};
use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::sync::Arc;

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
    model: Arc<pdbiox::ModelCif>,
    index: usize,
}

#[pymethods]
impl PyModelCategory {
    #[getter]
    fn name(&self) -> &str {
        self.inner().map_or("", pdbiox::ModelCategory::name)
    }

    #[getter]
    fn items(&self) -> Vec<String> {
        self.inner().map_or_else(Vec::new, |category| {
            category.items().iter().map(ToString::to_string).collect()
        })
    }

    #[getter]
    fn rows(&self) -> Vec<PyModelRow> {
        let Some(category) = self.inner() else {
            return Vec::new();
        };
        category
            .rows()
            .map(|row| PyModelRow {
                values: category
                    .items()
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        row.value(index).map_or_else(
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
        self.inner()?
            .value(item, row)
            .map(|inner| PyCifValue { inner })
    }
}

impl PyModelCategory {
    fn inner(&self) -> Option<&pdbiox::ModelCategory> {
        self.model.categories.get(self.index)
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
    model: Arc<pdbiox::ModelCif>,
}

#[pymethods]
impl PyQualityMetrics {
    #[getter]
    fn definitions(&self) -> Vec<PyMetricDefinition> {
        self.model
            .confidence()
            .definitions()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn global_metrics(&self) -> Vec<PyGlobalMetric> {
        self.model.confidence().global().map(Into::into).collect()
    }

    #[getter]
    fn local(&self) -> Vec<PyLocalMetric> {
        self.model.confidence().local().map(Into::into).collect()
    }

    #[getter]
    fn pairwise(&self) -> Vec<PyPairwiseMetric> {
        self.model.confidence().pairwise().map(Into::into).collect()
    }

    fn plddt(&self) -> Vec<PyLocalMetric> {
        self.model.confidence().plddt().map(Into::into).collect()
    }

    fn pae(&self) -> Vec<PyPairwiseMetric> {
        self.model.confidence().pae().map(Into::into).collect()
    }

    fn ptm(&self) -> Vec<PyGlobalMetric> {
        self.model.confidence().ptm().map(Into::into).collect()
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

metadata_projection!(PyModelDescription, "ModelDescription", pdbiox::ModelDescription<'_> {
    id: String,
    assembly_id: Option<String>,
    name: Option<String>,
    model_type: Option<String>,
});
metadata_projection!(PyTarget, "Target", pdbiox::Target<'_> {
    entity_id: String,
    data_id: Option<String>,
    origin: Option<String>,
});
metadata_projection!(PyTemplate, "Template", pdbiox::Template<'_> {
    id: String,
    target_chain_id: Option<String>,
    name: Option<String>,
    origin: Option<String>,
    entity_type: Option<String>,
});
metadata_projection!(PyProtocolStep, "ProtocolStep", pdbiox::ProtocolStep<'_> {
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
    pub(crate) inner: Arc<pdbiox::ModelCif>,
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
            .enumerate()
            .map(|(index, _)| PyModelCategory {
                model: Arc::clone(&self.inner),
                index,
            })
            .collect()
    }

    #[getter]
    fn models(&self) -> Vec<PyModelDescription> {
        self.inner.models().map(Into::into).collect()
    }

    #[getter]
    fn targets(&self) -> Vec<PyTarget> {
        self.inner.targets().map(Into::into).collect()
    }

    #[getter]
    fn templates(&self) -> Vec<PyTemplate> {
        self.inner.templates().map(Into::into).collect()
    }

    #[getter]
    fn protocol_steps(&self) -> Vec<PyProtocolStep> {
        self.inner.protocol_steps().map(Into::into).collect()
    }

    #[getter]
    fn software_groups(&self) -> Vec<PySoftwareGroup> {
        self.inner.software_groups().map(Into::into).collect()
    }

    #[getter]
    fn quality(&self) -> PyQualityMetrics {
        PyQualityMetrics {
            model: Arc::clone(&self.inner),
        }
    }

    fn confidence(&self) -> PyQualityMetrics {
        PyQualityMetrics::from_model(Arc::clone(&self.inner))
    }
}

#[pymethods]
impl PyStructure {
    fn model_cif(&self) -> Option<PyModelCif> {
        self.structure()
            .extensions()
            .get_shared::<pdbiox::ModelCif>(pdbiox::MODEL_CIF_EXTENSION)
            .map(|inner| PyModelCif { inner })
    }

    fn confidence(&self) -> Option<PyQualityMetrics> {
        self.model_cif()
            .map(|model| PyQualityMetrics::from_model(model.inner))
    }
}

#[pyfunction]
pub(crate) fn lower_model_cif(py: Python<'_>, document: &PyCifDocument) -> PyResult<PyModelCif> {
    // The lowering itself never touches the interpreter; only reporting does.
    let (model, findings) = py
        .detach(|| pdbiox::modelcif::lower(&document.inner))
        .map_err(|error| pyo3::exceptions::PyMemoryError::new_err(error.to_string()))?;
    if findings.is_empty() {
        Ok(PyModelCif {
            inner: Arc::new(model),
        })
    } else {
        Err(read_error(py, &findings))
    }
}

impl PyQualityMetrics {
    fn from_model(model: Arc<pdbiox::ModelCif>) -> Self {
        Self { model }
    }
}

impl From<pdbiox::MetricDefinition<'_>> for PyMetricDefinition {
    fn from(value: pdbiox::MetricDefinition<'_>) -> Self {
        Self {
            id: value.id.to_string(),
            name: value.name.map(ToString::to_string),
            metric_type: value.metric_type.to_string(),
            mode: value.mode.to_string(),
            software_group_id: value.software_group_id.map(ToString::to_string),
        }
    }
}

impl From<pdbiox::GlobalMetric<'_>> for PyGlobalMetric {
    fn from(value: pdbiox::GlobalMetric<'_>) -> Self {
        Self {
            model_id: value.model_id.to_string(),
            metric_id: value.metric_id.to_string(),
            value: value.value,
        }
    }
}

impl From<pdbiox::LocalMetric<'_>> for PyLocalMetric {
    fn from(value: pdbiox::LocalMetric<'_>) -> Self {
        Self {
            model_id: value.model_id.to_string(),
            chain_id: value.chain_id.to_string(),
            sequence_id: value.sequence_id,
            component_id: value.component_id.map(ToString::to_string),
            metric_id: value.metric_id.to_string(),
            value: value.value,
        }
    }
}

impl From<pdbiox::PairwiseMetric<'_>> for PyPairwiseMetric {
    fn from(value: pdbiox::PairwiseMetric<'_>) -> Self {
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

impl From<pdbiox::ModelDescription<'_>> for PyModelDescription {
    fn from(value: pdbiox::ModelDescription<'_>) -> Self {
        Self {
            id: value.id.to_string(),
            assembly_id: value.assembly_id.map(ToString::to_string),
            name: value.name.map(ToString::to_string),
            model_type: value.model_type.map(ToString::to_string),
        }
    }
}

impl From<pdbiox::Target<'_>> for PyTarget {
    fn from(value: pdbiox::Target<'_>) -> Self {
        Self {
            entity_id: value.entity_id.to_string(),
            data_id: value.data_id.map(ToString::to_string),
            origin: value.origin.map(ToString::to_string),
        }
    }
}

impl From<pdbiox::Template<'_>> for PyTemplate {
    fn from(value: pdbiox::Template<'_>) -> Self {
        Self {
            id: value.id.to_string(),
            target_chain_id: value.target_chain_id.map(ToString::to_string),
            name: value.name.map(ToString::to_string),
            origin: value.origin.map(ToString::to_string),
            entity_type: value.entity_type.map(ToString::to_string),
        }
    }
}

impl From<pdbiox::ProtocolStep<'_>> for PyProtocolStep {
    fn from(value: pdbiox::ProtocolStep<'_>) -> Self {
        Self {
            protocol_id: value.protocol_id.to_string(),
            step_id: value.step_id.to_string(),
            method_type: value.method_type.to_string(),
            name: value.name.map(ToString::to_string),
            details: value.details.map(ToString::to_string),
            software_group_id: value.software_group_id.map(ToString::to_string),
        }
    }
}

impl From<pdbiox::SoftwareGroup<'_>> for PySoftwareGroup {
    fn from(value: pdbiox::SoftwareGroup<'_>) -> Self {
        Self {
            group: value.group_id.to_string(),
            software: value.software_id.to_string(),
            parameter_group: value.parameter_group_id.map(ToString::to_string),
        }
    }
}
