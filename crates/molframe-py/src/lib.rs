//! Curated Python bindings for the `MolFrame` facade.
//!
//! Python mirrors the stable Rust contract rather than the workspace crate
//! graph. Native storage and kernels remain in Rust; this crate owns only
//! lifetime-safe Python views and small conversion boundaries.

#![deny(unsafe_op_in_unsafe_fn)]

mod analysis_result;
#[cfg(any(feature = "analysis", feature = "validation", feature = "spatial"))]
mod backend;
mod bindings;
mod catalog;
#[cfg(feature = "chemistry")]
mod chemistry;
#[cfg(feature = "compare")]
mod compare;
#[cfg(feature = "crystal")]
mod crystal;
#[cfg(feature = "crystal")]
mod crystal_statistics;
mod formats;
#[cfg(feature = "analysis")]
mod governed;
mod hierarchy;
mod native_source;
mod policy;
#[cfg(feature = "query")]
mod query_aliases;
mod query_cache;
mod query_messages;
#[cfg(feature = "query")]
mod selection_expr;
#[cfg(feature = "sequence")]
mod sequence;
#[cfg(feature = "spatial")]
mod spatial;
#[cfg(feature = "surface")]
mod surface;
mod table;
#[cfg(feature = "trajectory")]
mod trajectory;
#[cfg(feature = "validation")]
mod validation;
#[cfg(feature = "analysis")]
mod workflow;

#[cfg(feature = "analysis")]
use bindings::{PyContactTable, atom_contacts};
use bindings::{PyQuery, PyReader, PySelection, PyStructure, PyStructureEditor, read};
#[cfg(feature = "geometry")]
use bindings::{centroid, distance_matrix, rmsd};
use hierarchy::{
    PyAtom, PyAtoms, PyChain, PyChains, PyModel, PyModels, PyResidue, PyResidueSelection,
    PyResidues,
};
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
    // The version of the crate this binary was compiled from, which is the
    // version of the wheel: maturin takes the package version from it.
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    module.add_class::<PyStructure>()?;
    module.add_class::<PyStructureEditor>()?;
    module.add_class::<PyAtom>()?;
    module.add_class::<PyAtoms>()?;
    module.add_class::<PyResidue>()?;
    module.add_class::<PyResidueSelection>()?;
    module.add_class::<PyResidues>()?;
    module.add_class::<PyChain>()?;
    module.add_class::<PyChains>()?;
    module.add_class::<PyModel>()?;
    module.add_class::<PyModels>()?;
    module.add_class::<PySelection>()?;
    module.add_class::<policy::PyAnalysisPolicy>()?;
    analysis_result::register(module)?;
    module.add_class::<table::PyTable>()?;
    module.add_class::<PyQuery>()?;
    module.add(
        "QueryError",
        module.py().get_type::<query_messages::QueryError>(),
    )?;
    module.add(
        "QueryWarning",
        module.py().get_type::<query_messages::QueryWarning>(),
    )?;
    #[cfg(feature = "query")]
    module.add_class::<query_aliases::PyQueryAliases>()?;
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
        analysis.add_function(wrap_pyfunction!(bindings::contacts, &analysis)?)?;
        analysis.add_function(wrap_pyfunction!(governed::hydrogen_bonds, &analysis)?)?;
        analysis.add_function(wrap_pyfunction!(governed::salt_bridges, &analysis)?)?;
        analysis.add("ContactTable", module.getattr("ContactTable")?)?;
    }
    module.add_submodule(&analysis)?;

    let selection = PyModule::new(py, "sel")?;
    #[cfg(feature = "query")]
    selection_expr::register(&selection)?;
    module.add_submodule(&selection)?;

    let query = PyModule::new(py, "query")?;
    #[cfg(feature = "query")]
    query.add_function(wrap_pyfunction!(query_aliases::complete, &query)?)?;
    module.add_submodule(&query)?;

    let crystal = PyModule::new(py, "crystal")?;
    #[cfg(feature = "crystal")]
    crystal::register(&crystal)?;
    module.add_submodule(&crystal)?;

    let compare = PyModule::new(py, "compare")?;
    #[cfg(feature = "compare")]
    compare::register(&compare)?;
    module.add_submodule(&compare)?;

    let surface = PyModule::new(py, "surface")?;
    #[cfg(feature = "surface")]
    surface::register(&surface)?;
    module.add_submodule(&surface)?;

    let sequence = PyModule::new(py, "sequence")?;
    #[cfg(feature = "sequence")]
    sequence::register(&sequence)?;
    module.add_submodule(&sequence)?;

    let chemistry = PyModule::new(py, "chemistry")?;
    #[cfg(feature = "chemistry")]
    chemistry::register(&chemistry)?;
    module.add_submodule(&chemistry)?;

    let validation = PyModule::new(py, "validation")?;
    #[cfg(feature = "validation")]
    validation::register(&validation)?;
    module.add_submodule(&validation)?;

    let spatial = PyModule::new(py, "spatial")?;
    #[cfg(feature = "spatial")]
    spatial::register(&spatial)?;
    module.add_submodule(&spatial)?;

    let trajectory = PyModule::new(py, "trajectory")?;
    #[cfg(feature = "trajectory")]
    trajectory::register(&trajectory)?;
    module.add_submodule(&trajectory)?;

    module.add_submodule(&PyModule::new(py, "motif")?)?;
    let formats = PyModule::new(py, "formats")?;
    formats::register(&formats)?;
    module.add_submodule(&formats)
}

#[cfg(test)]
#[path = "module_tests.rs"]
mod tests;
