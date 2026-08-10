//! Migration entry points that delegate to the native Rust facade.

use crate::atom::PyAtoms;
use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::exceptions::PyUserWarning;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

static BIOPYTHON_NOTE_EMITTED: AtomicBool = AtomicBool::new(false);
static MDANALYSIS_NOTE_EMITTED: AtomicBool = AtomicBool::new(false);

type NativeReader = for<'py> fn(Python<'py>, PathBuf) -> PyResult<PyStructure>;

/// Biopython-shaped parser returning a native pdbiox structure.
#[pyclass(name = "PDBParser", module = "pdbiox.compat.biopython", frozen)]
pub(crate) struct PyPdbParser {
    reader: NativeReader,
}

#[pymethods]
impl PyPdbParser {
    #[new]
    fn new(py: Python<'_>) -> PyResult<Self> {
        warn_once(
            py,
            &BIOPYTHON_NOTE_EMITTED,
            c"pdbiox.compat.biopython.PDBParser delegates to pdbiox.read and returns pdbiox.Structure",
        )?;
        Ok(Self {
            reader: load_structure,
        })
    }

    /// Reads through the Rust format dispatcher; the compatibility identifier
    /// is accepted for call-site migration but does not rename stored data.
    fn get_structure(
        &self,
        py: Python<'_>,
        _identifier: &str,
        path: PathBuf,
    ) -> PyResult<PyStructure> {
        (self.reader)(py, path)
    }
}

/// MDAnalysis-shaped single-topology container over native pdbiox handles.
#[pyclass(name = "Universe", module = "pdbiox.compat.mdanalysis", frozen)]
pub(crate) struct PyUniverse {
    structure: PyStructure,
}

#[pymethods]
impl PyUniverse {
    #[new]
    fn new(py: Python<'_>, topology: PathBuf) -> PyResult<Self> {
        warn_once(
            py,
            &MDANALYSIS_NOTE_EMITTED,
            c"pdbiox.compat.mdanalysis.Universe delegates topology reading to pdbiox.read; use pdbiox.Structure directly",
        )?;
        load_structure(py, topology).map(|structure| Self { structure })
    }

    /// Native structure snapshot held by this compatibility container.
    #[getter]
    fn structure(&self) -> PyStructure {
        self.structure.clone()
    }

    /// Native atom collection from the same immutable structure snapshot.
    #[getter]
    fn atoms(&self) -> PyAtoms {
        PyAtoms::new(self.structure.structure().clone())
    }
}

fn load_structure(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    pdbiox::read(path)
        .map(PyStructure::new)
        .map_err(|findings| read_error(py, &findings))
}

fn warn_once(py: Python<'_>, emitted: &AtomicBool, message: &std::ffi::CStr) -> PyResult<()> {
    if emitted.swap(true, Ordering::Relaxed) {
        return Ok(());
    }
    PyErr::warn(py, &py.get_type::<PyUserWarning>(), message, 2)
}

#[cfg(test)]
#[path = "compatibility_tests.rs"]
mod tests;
