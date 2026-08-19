//! Owned adapters for CCD-derived side-chain paths.

use super::components::PyComponent;
use ::pdbiox;
use pyo3::prelude::*;

#[pyclass(name = "SideChainRoles", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySideChainRoles {
    #[pyo3(get)]
    nitrogen: String,
    #[pyo3(get)]
    alpha_carbon: String,
    #[pyo3(get)]
    side_chain_atoms: Vec<String>,
}

#[pymethods]
impl PySideChainRoles {
    #[new]
    fn new(nitrogen: String, alpha_carbon: String, side_chain_atoms: Vec<String>) -> Self {
        Self {
            nitrogen,
            alpha_carbon,
            side_chain_atoms,
        }
    }
}

#[pyclass(name = "SideChainDefinition", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySideChainDefinition(pub(crate) pdbiox::SideChainDefinition);

#[pymethods]
impl PySideChainDefinition {
    #[getter]
    fn atoms(&self) -> Vec<String> {
        self.0.atoms.iter().map(ToString::to_string).collect()
    }

    fn torsion_count(&self) -> usize {
        self.0.torsion_count()
    }
}

#[pyfunction]
pub(crate) fn side_chain_definition(
    component: &PyComponent,
    roles: &PySideChainRoles,
) -> Option<PySideChainDefinition> {
    let atoms = roles
        .side_chain_atoms
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let roles = pdbiox::SideChainRoles {
        nitrogen: &roles.nitrogen,
        alpha_carbon: &roles.alpha_carbon,
        side_chain_atoms: &atoms,
    };
    pdbiox::side_chain_definition(&component.0, &roles).map(PySideChainDefinition)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySideChainRoles>()?;
    module.add_class::<PySideChainDefinition>()?;
    module.add_function(wrap_pyfunction!(side_chain_definition, module)?)?;
    Ok(())
}
