//! Curated Python bindings for the `MolFrame` facade.
//!
//! Python mirrors the stable Rust contract rather than the workspace crate
//! graph. Native storage and kernels remain in Rust; this crate owns only
//! lifetime-safe Python views and small conversion boundaries.

#![deny(unsafe_op_in_unsafe_fn)]

mod bindings;
mod catalog;
mod hierarchy;
mod workflow;

use bindings::{
    PyContactTable, PyQuery, PyReader, PySelection, PyStructure, PyStructureEditor, atom_contacts,
    centroid, distance_matrix, read, rmsd,
};
use hierarchy::{PyAtom, PyAtoms, PyChain, PyChains, PyModel, PyModels, PyResidue, PyResidues};
use pyo3::prelude::*;
use workflow::{PyCompiledWorkflow, PyWorkflow, PyWorkflowNode};

#[pymodule]
#[pyo3(name = "_native")]
fn native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStructure>()?;
    module.add_class::<PyStructureEditor>()?;
    module.add_class::<PyAtom>()?;
    module.add_class::<PyAtoms>()?;
    module.add_class::<PyResidue>()?;
    module.add_class::<PyResidues>()?;
    module.add_class::<PyChain>()?;
    module.add_class::<PyChains>()?;
    module.add_class::<PyModel>()?;
    module.add_class::<PyModels>()?;
    module.add_class::<PySelection>()?;
    module.add_class::<PyQuery>()?;
    module.add_class::<PyReader>()?;
    module.add_class::<PyContactTable>()?;
    module.add_class::<PyWorkflow>()?;
    module.add_class::<PyWorkflowNode>()?;
    module.add_class::<PyCompiledWorkflow>()?;
    module.add_function(wrap_pyfunction!(read, module)?)?;
    register_namespaces(module)?;
    catalog::validate_registration(module)
}

fn register_namespaces(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let geometry = PyModule::new(py, "geometry")?;
    geometry.add_function(wrap_pyfunction!(centroid, &geometry)?)?;
    geometry.add_function(wrap_pyfunction!(distance_matrix, &geometry)?)?;
    geometry.add_function(wrap_pyfunction!(rmsd, &geometry)?)?;
    module.add_submodule(&geometry)?;

    let analysis = PyModule::new(py, "analysis")?;
    analysis.add_function(wrap_pyfunction!(atom_contacts, &analysis)?)?;
    analysis.add("ContactTable", module.getattr("ContactTable")?)?;
    module.add_submodule(&analysis)?;

    for name in [
        "trajectory",
        "sequence",
        "crystal",
        "validation",
        "motif",
        "chemistry",
        "compare",
        "query",
        "spatial",
        "surface",
    ] {
        module.add_submodule(&PyModule::new(py, name)?)?;
    }
    let formats = PyModule::new(py, "formats")?;
    for name in ["cif", "bcif", "pdb", "modelcif"] {
        formats.add_submodule(&PyModule::new(py, name)?)?;
    }
    module.add_submodule(&formats)
}

#[cfg(test)]
#[path = "module_tests.rs"]
mod tests;
