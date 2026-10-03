//! Conversion-only cross-instance covalent link adapter.

use crate::bindings::PyStructure;
use molframe::crystal::AssemblyExt as _;
use pyo3::{exceptions::PyValueError, prelude::*};

#[derive(Clone)]
#[pyclass(
    name = "AssemblyBond",
    module = "molframe.crystal",
    frozen,
    skip_from_py_object
)]
struct PyAssemblyBond {
    #[pyo3(get)]
    first_instance: u32,
    #[pyo3(get)]
    first_atom: u32,
    #[pyo3(get)]
    second_instance: u32,
    #[pyo3(get)]
    second_atom: u32,
    #[pyo3(get)]
    order: &'static str,
}

/// Geometrically perceived covalent links between distinct assembly instances.
#[pyfunction]
#[pyo3(signature = (structure, id, *, model=0))]
fn assembly_covalent_links(
    py: Python<'_>,
    structure: &PyStructure,
    id: String,
    model: u32,
) -> PyResult<Vec<PyAssemblyBond>> {
    let source = structure.inner.clone();
    let links = py.detach(move || {
        let view = source
            .assembly(&id)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        view.covalent_links(
            molframe::ModelIndex::new(model),
            &molframe::ExecutionContext::default(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))
    })?;
    Ok(links
        .into_iter()
        .map(|bond| PyAssemblyBond {
            first_instance: bond.first_instance.get(),
            first_atom: bond.first_atom.get(),
            second_instance: bond.second_instance.get(),
            second_atom: bond.second_atom.get(),
            order: "single",
        })
        .collect())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAssemblyBond>()?;
    module.add_function(wrap_pyfunction!(assembly_covalent_links, module)?)
}
