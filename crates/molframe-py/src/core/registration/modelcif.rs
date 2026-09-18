//! Registration for typed `ModelCIF` projections.

use crate::modelcif::{
    PyGlobalMetric, PyLocalMetric, PyMetricDefinition, PyModelCategory, PyModelCif,
    PyModelDescription, PyModelRow, PyPairwiseMetric, PyProtocolStep, PyQualityMetrics,
    PySoftwareGroup, PyTarget, PyTemplate, lower_model_cif,
};
use pyo3::prelude::*;

pub(super) fn register_modelcif(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyModelRow>()?;
    module.add_class::<PyModelCategory>()?;
    module.add_class::<PyMetricDefinition>()?;
    module.add_class::<PyGlobalMetric>()?;
    module.add_class::<PyLocalMetric>()?;
    module.add_class::<PyPairwiseMetric>()?;
    module.add_class::<PyQualityMetrics>()?;
    module.add_class::<PyModelDescription>()?;
    module.add_class::<PyTarget>()?;
    module.add_class::<PyTemplate>()?;
    module.add_class::<PyProtocolStep>()?;
    module.add_class::<PySoftwareGroup>()?;
    module.add_class::<PyModelCif>()?;
    module.add_function(wrap_pyfunction!(lower_model_cif, module)?)?;
    crate::modelcif_write::register(module)?;
    Ok(())
}
