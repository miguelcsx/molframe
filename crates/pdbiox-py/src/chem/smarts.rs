//! Native SMARTS parsing and substructure matching adapters.

use super::components::PyComponent;
use crate::structure::PyStructure;
use ::pdbiox;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, SmartsError, PyValueError);
create_exception!(_native, SmartsDataError, PyValueError);

#[pyclass(name = "SmartsMatch", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySmartsMatch {
    #[pyo3(get)]
    atom_indices: Vec<usize>,
}

#[pyclass(name = "SmartsPattern", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySmartsPattern(pub(crate) pdbiox::SmartsPattern);

#[pymethods]
impl PySmartsPattern {
    #[new]
    fn new(text: &str) -> PyResult<Self> {
        pdbiox::SmartsPattern::parse(text)
            .map(Self)
            .map_err(|error| SmartsError::new_err(error.to_string()))
    }

    #[staticmethod]
    fn parse(text: &str) -> PyResult<Self> {
        Self::new(text)
    }

    fn find_matches(&self, py: Python<'_>, component: &PyComponent) -> Vec<PySmartsMatch> {
        let pattern = self.0.clone();
        let component = component.0.clone();
        py.detach(move || pattern.find_matches(&component))
            .into_iter()
            .map(|value| PySmartsMatch {
                atom_indices: value.atom_indices.into_vec(),
            })
            .collect()
    }

    fn find_structure_matches(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
    ) -> PyResult<Vec<PySmartsMatch>> {
        let pattern = self.0.clone();
        let structure = structure.structure().clone();
        py.detach(move || pattern.find_structure_matches(&structure))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| PySmartsMatch {
                        atom_indices: value.atom_indices.into_vec(),
                    })
                    .collect()
            })
            .map_err(|error| SmartsDataError::new_err(error.to_string()))
    }

    fn matches(&self, py: Python<'_>, component: &PyComponent) -> bool {
        let pattern = self.0.clone();
        let component = component.0.clone();
        py.detach(move || pattern.matches(&component))
    }

    fn atom_count(&self) -> usize {
        self.0.atom_count()
    }
}

#[pyfunction]
pub(crate) fn parse_smarts(py: Python<'_>, text: &str) -> PyResult<PySmartsPattern> {
    py.detach(move || -> PyResult<PySmartsPattern> { PySmartsPattern::new(text) })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("SmartsError", module.py().get_type::<SmartsError>())?;
    module.add("SmartsDataError", module.py().get_type::<SmartsDataError>())?;
    module.add_class::<PySmartsMatch>()?;
    module.add_class::<PySmartsPattern>()?;
    module.add_function(wrap_pyfunction!(parse_smarts, module)?)?;
    Ok(())
}
