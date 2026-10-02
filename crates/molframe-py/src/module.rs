//! The native module and the namespaces it registers.

use super::{
    analysis_result, bindings, catalog, chemistry, compare, crystal, formats, governed, hierarchy,
    policy, query_aliases, query_messages, selection_expr, sequence, spatial, surface, table,
    trajectory, validation, workflow,
};
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

#[pymodule]
#[pyo3(name = "_native")]
pub(super) fn native(module: &Bound<'_, PyModule>) -> PyResult<()> {
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

pub(super) fn register_namespaces(module: &Bound<'_, PyModule>) -> PyResult<()> {
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

    let motif = PyModule::new(py, "motif")?;
    #[cfg(feature = "motif")]
    crate::motif::register(&motif)?;
    module.add_submodule(&motif)?;
    let formats = PyModule::new(py, "formats")?;
    formats::register(&formats)?;
    module.add_submodule(&formats)
}
