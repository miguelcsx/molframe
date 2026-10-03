//! Conversion-only structure charge adapter.

use crate::bindings::{PyStructure, findings_error};
use molframe::chemistry as chem;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::path::PathBuf;

#[pyclass(
    name = "PartialCharges",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PyPartialCharges {
    #[pyo3(get)]
    values: Vec<f64>,
    #[pyo3(get)]
    source: &'static str,
    #[pyo3(get)]
    dictionary_version: Option<String>,
    #[pyo3(get)]
    parameter_profile: Option<&'static str>,
}

impl From<chem::PartialCharges> for PyPartialCharges {
    fn from(charges: chem::PartialCharges) -> Self {
        let (source, dictionary_version, parameter_profile) = match charges.source {
            chem::ChargeSource::File => ("file", None, None),
            chem::ChargeSource::Peoe {
                dictionary,
                options,
            } => (
                "peoe",
                Some(dictionary.as_str().to_owned()),
                Some(options.profile.name()),
            ),
        };
        Self {
            values: charges.values,
            source,
            dictionary_version,
            parameter_profile,
        }
    }
}

/// Stable atom-order partial charges, preferring complete file charge columns.
#[pyfunction]
#[pyo3(signature = (structure, components=None, *, version="unversioned"))]
fn partial_charges(
    py: Python<'_>,
    structure: &PyStructure,
    components: Option<PathBuf>,
    version: &str,
) -> PyResult<PyPartialCharges> {
    let source = structure.inner.clone();
    let version = molframe::DictionaryVersion::new(version);
    let result = py.detach(move || {
        let calculate = |provider: &dyn chem::ComponentProvider| {
            chem::partial_charges(source.engine(), provider, chem::PeoeOptions::default())
                .map_err(|error| PyValueError::new_err(error.to_string()))
        };
        if let Some(path) = components {
            let provider = molframe::read_component_dictionary(&path, version)
                .map_err(|findings| findings_error(&findings))?
                .0;
            calculate(&provider)
        } else {
            let provider = chem::MemoryProvider::new(version, [])
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
            calculate(&provider)
        }
    })?;
    Ok(result.into())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPartialCharges>()?;
    module.add_function(wrap_pyfunction!(partial_charges, module)?)
}
