//! Explicit biological-assembly structures, whole or in bounded pieces.

use crate::bindings::PyStructure;
use molframe::crystal::{AssemblyChunks, AssemblyExt as _, AssemblyView};
use pyo3::prelude::*;
use std::sync::Mutex;

fn view_of(structure: &PyStructure, id: &str) -> PyResult<AssemblyView> {
    let engine = structure.inner.engine();
    let set = engine
        .assembly_set()
        .ok_or_else(|| crate::error::value("the structure declares no biological assemblies"))?;
    AssemblyView::new(engine, set, id).map_err(crate::error::kernel)
}

/// The biological assembly `id` as one structure with every copy placed.
///
/// Copies get unique chain labels (`A`, `A-2`, ...) and an instance id
/// annotation. For capsid-scale assemblies use `assembly_chunks`, which holds
/// one bounded piece at a time.
#[pyfunction]
fn assembly_structure(py: Python<'_>, structure: &PyStructure, id: &str) -> PyResult<PyStructure> {
    let view = view_of(structure, id)?;
    let built = py
        .detach(|| view.materialize())
        .map_err(|findings| crate::error::from_findings(&findings, "assembly expansion failed"))?;
    Ok(PyStructure::new(molframe::Structure::from(built)))
}

/// Iterator over an assembly's pieces, each a structure of whole placements.
#[pyclass(name = "AssemblyChunks", module = "molframe", skip_from_py_object)]
struct PyAssemblyChunks {
    inner: Mutex<AssemblyChunks>,
}

#[pymethods]
impl PyAssemblyChunks {
    const fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyStructure>> {
        let next = py.detach(|| {
            self.inner
                .lock()
                .map_err(|_| crate::error::internal("assembly iterator lock poisoned"))
                .map(|mut chunks| chunks.next())
        })?;
        match next {
            None => Ok(None),
            Some(Ok(piece)) => Ok(Some(PyStructure::new(molframe::Structure::from(piece)))),
            Some(Err(findings)) => Err(crate::error::from_findings(
                &findings,
                "assembly expansion failed",
            )),
        }
    }
}

/// Stream assembly `id` as structures of about `max_atoms` atoms each.
///
/// A placement is never split, so a piece may exceed the bound by at most one
/// placement. Labels and atom ids continue across pieces exactly as in
/// `assembly_structure`.
#[pyfunction]
#[pyo3(signature = (structure, id, max_atoms = 1_000_000))]
fn assembly_chunks(
    structure: &PyStructure,
    id: &str,
    max_atoms: usize,
) -> PyResult<PyAssemblyChunks> {
    let view = view_of(structure, id)?;
    Ok(PyAssemblyChunks {
        inner: Mutex::new(view.chunks(max_atoms)),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAssemblyChunks>()?;
    module.add_function(wrap_pyfunction!(assembly_structure, module)?)?;
    module.add_function(wrap_pyfunction!(assembly_chunks, module)?)
}
