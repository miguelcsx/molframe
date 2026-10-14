//! Motion candidates preserve the caller's original topology and atom addresses.

use crate::bindings::PyStructure;
use crate::bindings::array::Array2;
use crate::execution::PyExecutionContext;
use crate::table::{PyTable, TableBuilder};
use molframe::analysis::{HydrogenBondCandidateOptions, hydrogen_bond_candidates_over_segment};
use numpy::PyArrayMethods;
use pyo3::prelude::*;

#[pyclass(
    name = "HydrogenBondCandidateOptions",
    module = "molframe.analysis",
    frozen,
    skip_from_py_object
)]
pub(crate) struct PyCandidateOptions {
    inner: HydrogenBondCandidateOptions,
}

#[pymethods]
impl PyCandidateOptions {
    #[new]
    #[pyo3(signature = (*, max_distance=3.5, maximum_candidates=2_000_000, backend="auto"))]
    fn new(max_distance: f32, maximum_candidates: usize, backend: &str) -> PyResult<Self> {
        let inner = HydrogenBondCandidateOptions {
            maximum_donor_acceptor_distance: max_distance,
            maximum_candidates,
            backend: crate::backend::parse(backend)?,
        };
        inner.validate().map_err(crate::error::kernel)?;
        Ok(Self { inner })
    }
}

/// Uses two borrowed C-contiguous float32 coordinate columns and existing CCD
/// roles. Atom addresses stay in the original source: no conformers are removed
/// or reindexed. Candidate triples require instantaneous geometric filtering.
#[pyfunction]
#[pyo3(signature = (structure, start, end, *, options=None, context=None))]
fn hydrogen_bond_candidates(
    py: Python<'_>,
    structure: &PyStructure,
    start: Array2<'_, '_, f32>,
    end: Array2<'_, '_, f32>,
    options: Option<&PyCandidateOptions>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyTable> {
    let source = structure.inner.clone();
    let start_view = start.readonly();
    let end_view = end.readonly();
    let start = crate::bindings::coordinates(&start_view)?;
    let end = crate::bindings::coordinates(&end_view)?;
    let options = match options {
        Some(options) => options.inner,
        None => PyCandidateOptions::new(3.5, 2_000_000, "auto")?.inner,
    };
    let candidates = crate::execution::run(py, context, |context| {
        hydrogen_bond_candidates_over_segment(source.engine(), start, end, options, context)
    })?
    .map_err(crate::error::kernel)?;
    let indices = |column: &[molframe::AtomIndex]| -> Vec<u32> {
        column.iter().map(|index| index.get()).collect()
    };
    Ok(TableBuilder::new(py, candidates.len())
        .indices("donor", &indices(candidates.donor()))
        .indices("hydrogen", &indices(candidates.hydrogen()))
        .indices("acceptor", &indices(candidates.acceptor()))
        .finish())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCandidateOptions>()?;
    module.add_function(wrap_pyfunction!(hydrogen_bond_candidates, module)?)
}
