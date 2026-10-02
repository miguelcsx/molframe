//! Mechanical adapter for Niggli cell reduction.

use molframe::UnitCell;
use molframe::crystal::{ReducedCell, niggli_reduce};
use pyo3::{exceptions::PyValueError, prelude::*};

/// A cell in its Niggli setting and the change of basis that reaches it.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ReducedCell",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyReducedCell(ReducedCell);

#[pymethods]
impl PyReducedCell {
    /// Edge lengths of the reduced cell in ångström.
    #[getter]
    const fn lengths(&self) -> [f64; 3] {
        self.0.lengths
    }

    /// Angles `α, β, γ` of the reduced cell in degrees.
    #[getter]
    const fn angles(&self) -> [f64; 3] {
        self.0.angles
    }

    /// Columns are the reduced edges written in the original edges.
    #[getter]
    const fn change_of_basis(&self) -> [[i32; 3]; 3] {
        self.0.change_of_basis
    }

    /// How many normalize-and-reduce rounds ran.
    #[getter]
    const fn iterations(&self) -> usize {
        self.0.iterations
    }

    /// Whether the last round found nothing left to reduce.
    #[getter]
    const fn converged(&self) -> bool {
        self.0.converged
    }

    fn __repr__(&self) -> String {
        format!(
            "ReducedCell(lengths={:?}, angles={:?})",
            self.0.lengths, self.0.angles
        )
    }
}

/// Reduces a primitive cell to its Niggli setting.
#[pyfunction]
#[pyo3(signature = (lengths, angles, *, epsilon=1e-9, iteration_limit=100))]
fn reduce_cell(
    lengths: [f64; 3],
    angles: [f64; 3],
    epsilon: f64,
    iteration_limit: usize,
) -> PyResult<PyReducedCell> {
    niggli_reduce(&UnitCell { lengths, angles }, epsilon, iteration_limit)
        .map(PyReducedCell)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyReducedCell>()?;
    module.add_function(wrap_pyfunction!(reduce_cell, module)?)
}
