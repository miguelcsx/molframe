//! Owned projection of preserved legacy PDB metadata.

use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "PdbHeaderRecord", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPdbHeaderRecord {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    line: String,
}

#[pyclass(name = "PdbHeaders", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPdbHeaders {
    records: Vec<PyPdbHeaderRecord>,
}

#[pymethods]
impl PyPdbHeaders {
    fn __len__(&self) -> usize {
        self.records.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    #[getter]
    fn records(&self) -> Vec<PyPdbHeaderRecord> {
        self.records.clone()
    }

    fn named(&self, name: &str) -> Vec<PyPdbHeaderRecord> {
        self.records
            .iter()
            .filter(|record| record.name == name)
            .cloned()
            .collect()
    }

    #[getter]
    fn classification(&self) -> Option<String> {
        self.records
            .iter()
            .find(|record| record.name == "HEADER")
            .and_then(|record| record.line.get(10..50))
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_owned)
    }

    #[getter]
    fn deposition_date(&self) -> Option<String> {
        self.records
            .iter()
            .find(|record| record.name == "HEADER")
            .and_then(|record| record.line.get(50..59))
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_owned)
    }
}

#[pymethods]
impl PyStructure {
    fn pdb_headers(&self) -> Option<PyPdbHeaders> {
        molframe::PdbHeadersExt::pdb_headers(self.structure()).map(|headers| PyPdbHeaders {
            records: headers
                .records()
                .iter()
                .map(|record| PyPdbHeaderRecord {
                    name: record.name().to_owned(),
                    line: record.line().to_owned(),
                })
                .collect(),
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPdbHeaderRecord>()?;
    module.add_class::<PyPdbHeaders>()?;
    Ok(())
}
