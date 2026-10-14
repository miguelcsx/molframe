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
///
/// Bonded (1-2), 1-3 and, unless `exclude_one_four=False`, 1-4 pairs are never
/// reported. With `perceive_bonds=True` the structure's bond table is completed
/// with intra-residue bonds, peptide and phosphodiester links and disulfides
/// first, so a structure read from a coordinate file is not flooded with its own
/// covalent neighbours.
#[pyfunction]
#[pyo3(signature = (structure, *, tolerance=0.4, radii=None, backend="auto", exclude_one_four=true, perceive_bonds=true, policy=None, context=None))]
fn clashes(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radii: Option<&str>,
    backend: &str,
    exclude_one_four: bool,
    perceive_bonds: bool,
    policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<crate::analysis_result::PyAnalysis> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(crate::error::value(
            "tolerance must be finite and non-negative",
        ));
    }
    let backend = crate::backend::parse(backend)?;
    let mut options = molframe::validation::ClashOptions::resolve(
        tolerance,
        &molframe::chemistry::TableRequest::<molframe::chemistry::VdwRadii>::choose(radii),
        backend,
        &molframe::chemistry::TableStore::global(),
    )
    .map_err(crate::error::kernel)?;
    options.exclude_one_four = exclude_one_four;
    options.perceive_bonds = perceive_bonds;
    let kernel = molframe::validation::clashes_with_kernel(options);
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
    crate::validation_builtin::register(module)?;
    module.add_class::<PyClashTable>()?;
    module.add_function(wrap_pyfunction!(clashes, module)?)
}
