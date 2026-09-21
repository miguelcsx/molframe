//! Curated Python bindings for the `MolFrame` facade.
//!
//! Python mirrors the stable Rust contract rather than the workspace crate
//! graph. Native storage and kernels remain in Rust; this crate owns only
//! lifetime-safe Python views and small conversion boundaries.

#![deny(unsafe_op_in_unsafe_fn)]

mod bindings;
mod catalog;
mod hierarchy;
mod native_source;
#[cfg(feature = "query")]
mod selection_expr;
#[cfg(feature = "analysis")]
mod workflow;

#[cfg(feature = "analysis")]
use bindings::{PyContactTable, atom_contacts};
use bindings::{PyQuery, PyReader, PySelection, PyStructure, PyStructureEditor, read};
#[cfg(feature = "geometry")]
use bindings::{centroid, distance_matrix, rmsd};
use hierarchy::{PyAtom, PyAtoms, PyChain, PyChains, PyModel, PyModels, PyResidue, PyResidues};
use pyo3::prelude::*;
#[cfg(feature = "analysis")]
use workflow::{PyCompiledWorkflow, PyWorkflow, PyWorkflowNode};

/// Borrows the native snapshot retained by a Python `molframe.Structure`.
///
/// Cloning the returned value only increments its shared storage ownership;
/// coordinate columns are not copied.
///
/// # Errors
///
/// Returns Python's type error when `object` is not a native `MolFrame`
/// structure.
pub fn structure_from_python(object: &Bound<'_, PyAny>) -> PyResult<molframe::Structure> {
    Ok(object
        .extract::<PyRef<'_, PyStructure>>()
        .map(|structure| structure.inner.clone())?)
}

pub use native_source::{NativeAtom, NativeBond, NativeStructureSource, NativeTopology};

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
    #[cfg(feature = "analysis")]
    module.add_class::<PyContactTable>()?;
    #[cfg(feature = "analysis")]
    module.add_class::<PyWorkflow>()?;
    #[cfg(feature = "analysis")]
    module.add_class::<PyWorkflowNode>()?;
    #[cfg(feature = "analysis")]
    module.add_class::<PyCompiledWorkflow>()?;
    module.add_function(wrap_pyfunction!(read, module)?)?;
    register_namespaces(module)?;
    catalog::validate_registration(module)
}

fn register_namespaces(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let geometry = PyModule::new(py, "geometry")?;
    #[cfg(feature = "geometry")]
    {
        geometry.add_function(wrap_pyfunction!(centroid, &geometry)?)?;
        geometry.add_function(wrap_pyfunction!(distance_matrix, &geometry)?)?;
        geometry.add_function(wrap_pyfunction!(rmsd, &geometry)?)?;
    }
    module.add_submodule(&geometry)?;

    let analysis = PyModule::new(py, "analysis")?;
    #[cfg(feature = "analysis")]
    {
        analysis.add_function(wrap_pyfunction!(atom_contacts, &analysis)?)?;
        analysis.add("ContactTable", module.getattr("ContactTable")?)?;
    }
    module.add_submodule(&analysis)?;

    let selection = PyModule::new(py, "sel")?;
    #[cfg(feature = "query")]
    selection_expr::register(&selection)?;
    module.add_submodule(&selection)?;

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
