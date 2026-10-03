//! Conversion-only affine potential-grid adapters.

use crate::bindings::PyStructure;
use molframe::analysis;
use pyo3::prelude::*;

#[pyclass(
    name = "GridSpec",
    module = "molframe.analysis",
    frozen,
    skip_from_py_object
)]
#[derive(Clone, Copy)]
struct PyGridSpec {
    inner: analysis::GridSpec,
}

#[pymethods]
impl PyGridSpec {
    #[new]
    fn new(voxel_to_world: [[f64; 4]; 4], dimensions: [usize; 3]) -> Self {
        Self {
            inner: analysis::GridSpec {
                voxel_to_world,
                dimensions,
            },
        }
    }

    #[getter]
    fn voxel_to_world(&self) -> [[f64; 4]; 4] {
        self.inner.voxel_to_world
    }

    #[getter]
    fn dimensions(&self) -> [usize; 3] {
        self.inner.dimensions
    }
}

#[pyclass(
    name = "ScalarGrid",
    module = "molframe.analysis",
    frozen,
    skip_from_py_object
)]
struct PyScalarGrid {
    #[pyo3(get)]
    spec: PyGridSpec,
    #[pyo3(get)]
    values: Vec<f64>,
}

/// Screened Coulomb contact field, in kT/e at 298 K, not Poisson–Boltzmann.
#[pyfunction]
#[pyo3(signature = (structure, charges, spec, *, cutoff=analysis::CONTACT_POTENTIAL_CUTOFF))]
fn contact_potential(
    py: Python<'_>,
    structure: &PyStructure,
    charges: Vec<f64>,
    spec: &PyGridSpec,
    cutoff: f32,
) -> PyResult<PyScalarGrid> {
    let source = structure.inner.clone();
    let spec = spec.inner;
    let field = py
        .detach(move || analysis::contact_potential(source.engine(), &charges, spec, cutoff))
        .map_err(crate::error::failure)?;
    Ok(PyScalarGrid {
        spec: PyGridSpec { inner: field.spec },
        values: field.values,
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGridSpec>()?;
    module.add_class::<PyScalarGrid>()?;
    module.add_function(wrap_pyfunction!(contact_potential, module)?)
}
