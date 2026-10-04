//! Mechanical adapters for structure validation reports.

use crate::bindings::PyStructure;
use numpy::{PyArray1, ToPyArray};
use pyo3::prelude::*;

/// Atom pairs whose van der Waals spheres interpenetrate.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ClashTable",
    frozen,
    skip_from_py_object,
    module = "molframe.validation"
)]
struct PyClashTable {
    table: molframe::validation::ClashTable,
}

#[pymethods]
impl PyClashTable {
    fn __len__(&self) -> usize {
        self.table.len()
    }

    /// Lower atom index of each pair.
    #[getter]
    fn first<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.table
            .first()
            .iter()
            .map(|index| index.get())
            .collect::<Vec<_>>()
            .to_pyarray(py)
    }

    /// Higher atom index of each pair.
    #[getter]
    fn second<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.table
            .second()
            .iter()
            .map(|index| index.get())
            .collect::<Vec<_>>()
            .to_pyarray(py)
    }

    /// Interpenetration depth of each pair in ångström.
    #[getter]
    fn overlap<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        self.table.overlaps().to_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!("ClashTable(pairs={})", self.table.len())
    }
}

/// Steric clashes under an explicit policy: pairs overlapping by more than
/// `tolerance` ångström, with status, coverage and provenance.
#[pyfunction]
#[pyo3(signature = (structure, *, tolerance=0.4, radii="bondi", backend="auto", policy=None, context=None))]
fn clashes(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radii: &str,
    backend: &str,
    policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<crate::analysis_result::PyAnalysis> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(crate::error::value(
            "tolerance must be finite and non-negative",
        ));
    }
    let set = crate::chemistry::radius_set(radii)?;
    let backend = crate::backend::parse(backend)?;
    let kernel = molframe::validation::clashes_kernel(tolerance, set, backend);
    crate::governed::run(
        py,
        structure,
        &crate::policy::policy_of(policy),
        &kernel,
        |py, table| Ok(Py::new(py, PyClashTable { table })?.into_any()),
        context,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::validation_checks::register(module)?;
    module.add_class::<PyClashTable>()?;
    module.add_function(wrap_pyfunction!(clashes, module)?)
}
