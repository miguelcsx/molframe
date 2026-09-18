//! Native MOL, SDF and MOL2 parser/writer functions.

use super::molecule2::PyMol2Record;
use super::molecules::{Mol2Error, MolError, PyMolRecord};
use ::molframe::chem as molframe;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn parse_mol_record(py: Python<'_>, text: String) -> PyResult<PyMolRecord> {
    py.detach(move || molframe::parse_mol_record(&text))
        .map(super::molecules::PyMolRecord)
        .map_err(|error| MolError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_sdf_records(py: Python<'_>, text: String) -> PyResult<Vec<PyMolRecord>> {
    py.detach(move || molframe::parse_sdf_records(&text))
        .map(|records| {
            records
                .into_iter()
                .map(super::molecules::PyMolRecord)
                .collect()
        })
        .map_err(|error| MolError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_mol(py: Python<'_>, record: &PyMolRecord) -> PyResult<String> {
    let record = record.0.clone();
    py.detach(move || molframe::write_mol(&record))
        .map_err(|error| MolError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_sdf(py: Python<'_>, records: Vec<PyMolRecord>) -> PyResult<String> {
    let records = records
        .into_iter()
        .map(|record| record.0)
        .collect::<Vec<_>>();
    py.detach(move || molframe::write_sdf(&records))
        .map_err(|error| MolError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_mol2_record(py: Python<'_>, text: String) -> PyResult<PyMol2Record> {
    py.detach(move || molframe::parse_mol2_record(&text))
        .map(PyMol2Record)
        .map_err(|error| Mol2Error::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_mol2(py: Python<'_>, record: &PyMol2Record) -> PyResult<String> {
    let record = record.0.clone();
    py.detach(move || molframe::write_mol2(&record))
        .map_err(|error| Mol2Error::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(parse_mol_record, module)?)?;
    module.add_function(wrap_pyfunction!(parse_sdf_records, module)?)?;
    module.add_function(wrap_pyfunction!(write_mol, module)?)?;
    module.add_function(wrap_pyfunction!(write_sdf, module)?)?;
    module.add_function(wrap_pyfunction!(parse_mol2_record, module)?)?;
    module.add_function(wrap_pyfunction!(write_mol2, module)?)?;
    Ok(())
}
