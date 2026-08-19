//! Lossless fixed-column and Hybrid-36 primitives from the PDB facade.

use pyo3::prelude::*;

struct SelectionFilter {
    selection: Option<pdbiox::AtomSelection>,
}

impl pdbiox::Select for SelectionFilter {
    fn accept_atom(&self, atom: pdbiox::AtomIndex) -> bool {
        self.selection
            .as_ref()
            .is_none_or(|selection| selection.contains(atom.get()))
    }
}

#[pyfunction]
pub(crate) fn pdb_field_text(line: &str, start: usize, stop: usize) -> String {
    pdbiox::pdb::fixed::text(line, start, stop).to_owned()
}

#[pyfunction]
pub(crate) fn pdb_field_raw(line: &str, start: usize, stop: usize) -> String {
    pdbiox::pdb::fixed::raw(line, start, stop).to_owned()
}

#[pyfunction]
pub(crate) fn pdb_field_record(line: &str) -> String {
    pdbiox::pdb::fixed::record(line).to_owned()
}

#[pyfunction]
pub(crate) fn pdb_field_integer(line: &str, start: usize, stop: usize) -> Option<i64> {
    pdbiox::pdb::fixed::integer(line, start, stop)
}

#[pyfunction]
pub(crate) fn pdb_field_real(line: &str, start: usize, stop: usize) -> Option<f64> {
    pdbiox::pdb::fixed::real(line, start, stop)
}

#[pyfunction]
pub(crate) fn hybrid36_decode(field: &str, width: u32) -> Option<i64> {
    pdbiox::pdb::hybrid36::decode(field, width)
}

#[pyfunction]
pub(crate) fn hybrid36_encode(value: i64, width: u32) -> Option<String> {
    pdbiox::pdb::hybrid36::encode(value, width)
}

#[pyfunction]
pub(crate) fn hybrid36_needs_encoding(value: i64, width: u32) -> bool {
    pdbiox::pdb::hybrid36::needs_encoding(value, width)
}

#[pyfunction(name = "pdb_write_selected")]
#[pyo3(signature = (structure, options, selection=None))]
pub(crate) fn write_selected(
    py: Python<'_>,
    structure: &crate::structure::PyStructure,
    options: &crate::io::PyPdbWriteOptions,
    selection: Option<&crate::query::PySelection>,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let options = options.0.clone();
    let filter = SelectionFilter {
        selection: selection.map(|value| value.inner.clone()),
    };
    py.detach(move || pdbiox::pdb::write_selected(&structure, &options, &filter))
        .map_err(|findings| crate::errors::read_error(py, &findings))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("PDB_HEADERS_EXTENSION", pdbiox::pdb::PDB_HEADERS_EXTENSION)?;
    module.add(
        "MMTF_METADATA_EXTENSION",
        pdbiox::pdb::MMTF_METADATA_EXTENSION,
    )?;
    module.add_function(wrap_pyfunction!(pdb_field_text, module)?)?;
    module.add_function(wrap_pyfunction!(pdb_field_raw, module)?)?;
    module.add_function(wrap_pyfunction!(pdb_field_record, module)?)?;
    module.add_function(wrap_pyfunction!(pdb_field_integer, module)?)?;
    module.add_function(wrap_pyfunction!(pdb_field_real, module)?)?;
    module.add_function(wrap_pyfunction!(hybrid36_decode, module)?)?;
    module.add_function(wrap_pyfunction!(hybrid36_encode, module)?)?;
    module.add_function(wrap_pyfunction!(hybrid36_needs_encoding, module)?)?;
    module.add_function(wrap_pyfunction!(write_selected, module)?)?;
    Ok(())
}
