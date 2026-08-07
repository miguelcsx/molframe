//! Native-module registration and facade entry points.

use crate::atom::{PyAtom, PyAtoms};
use crate::bonds::PyBonds;
use crate::edit::PyCoordinateEdit;
use crate::errors::{read_error, register};
use crate::hierarchy::{
    PyChain, PyChains, PyModel, PyModels, PyResidue, PyResidueAtoms, PyResidues,
};
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyfunction]
fn read(py: Python<'_>, path: PathBuf) -> PyResult<PyStructure> {
    pdbiox::read(path)
        .map(PyStructure::new)
        .map_err(|findings| read_error(py, &findings))
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register(module)?;
    module.add_class::<PyStructure>()?;
    module.add_class::<PyCoordinateEdit>()?;
    module.add_class::<PyAtoms>()?;
    module.add_class::<PyAtom>()?;
    module.add_class::<PyBonds>()?;
    module.add_class::<PyModels>()?;
    module.add_class::<PyModel>()?;
    module.add_class::<PyChains>()?;
    module.add_class::<PyChain>()?;
    module.add_class::<PyResidues>()?;
    module.add_class::<PyResidue>()?;
    module.add_class::<PyResidueAtoms>()?;
    module.add_function(wrap_pyfunction!(read, module)?)?;
    Ok(())
}
