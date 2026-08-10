//! Nucleic torsion and sugar-pucker projections.

use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "NucleicTorsions", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNucleicTorsions {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    alpha: Option<f64>,
    #[pyo3(get)]
    beta: Option<f64>,
    #[pyo3(get)]
    gamma: Option<f64>,
    #[pyo3(get)]
    delta: Option<f64>,
    #[pyo3(get)]
    epsilon: Option<f64>,
    #[pyo3(get)]
    zeta: Option<f64>,
    #[pyo3(get)]
    chi: Option<f64>,
}

#[pyclass(name = "Pucker", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPucker {
    #[pyo3(get)]
    phase_degrees: f64,
    #[pyo3(get)]
    amplitude: f64,
}

#[pymethods]
impl PyStructure {
    fn nucleic_torsions(&self, py: Python<'_>) -> PyResult<Vec<PyNucleicTorsions>> {
        let structure = self.structure().clone();
        py.detach(move || pdbiox::analysis::nucleic_torsions(&structure))
            .map(|values| values.into_iter().map(PyNucleicTorsions::from).collect())
            .map_err(value_error)
    }
}

#[pyfunction]
pub(crate) fn sugar_pucker(py: Python<'_>, torsions: [f64; 5]) -> PyPucker {
    py.detach(move || pdbiox::analysis::sugar_pucker(torsions))
        .into()
}

impl From<pdbiox::analysis::NucleicTorsions> for PyNucleicTorsions {
    fn from(value: pdbiox::analysis::NucleicTorsions) -> Self {
        Self {
            residue: value.residue.get(),
            alpha: value.alpha,
            beta: value.beta,
            gamma: value.gamma,
            delta: value.delta,
            epsilon: value.epsilon,
            zeta: value.zeta,
            chi: value.chi,
        }
    }
}

impl From<pdbiox::analysis::Pucker> for PyPucker {
    fn from(value: pdbiox::analysis::Pucker) -> Self {
        Self {
            phase_degrees: value.phase_degrees,
            amplitude: value.amplitude,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
