//! Mechanical adapters for batched TM-score and its opt-in prefilter.

use crate::bindings::array::Readonly2;
use numpy::{IntoPyArray, PyArray1, PyUntypedArrayMethods};
use pyo3::prelude::*;

use crate::execution::PyExecutionContext;

/// Scores and the number of pairs each screen stage settled.
#[pyclass(name = "TmScreen", frozen, module = "molframe.compare")]
pub(crate) struct PyTmScreen {
    scores: Py<PyArray1<f64>>,
    evaluated: u64,
    pruned_bound: u64,
    pruned_coarse: u64,
}

#[pymethods]
impl PyTmScreen {
    /// One TM-score per requested pair; `nan` where the screen dropped it.
    #[getter]
    fn scores<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.scores.bind(py).clone()
    }

    /// Pairs that went through the full search.
    #[getter]
    fn evaluated(&self) -> u64 {
        self.evaluated
    }

    /// Pairs dropped by the exact upper bound.
    #[getter]
    fn pruned_bound(&self) -> u64 {
        self.pruned_bound
    }

    /// Pairs dropped by the heuristic coarse stage.
    #[getter]
    fn pruned_coarse(&self) -> u64 {
        self.pruned_coarse
    }
}

fn views<'a>(arrays: &'a [Readonly2<'_, f32>]) -> PyResult<Vec<&'a [[f32; 3]]>> {
    arrays
        .iter()
        .map(|array| crate::bindings::coordinates(array))
        .collect()
}

fn pair_list(pairs: &numpy::PyReadonlyArray2<'_, u64>) -> PyResult<Vec<(usize, usize)>> {
    if pairs.shape()[1] != 2 {
        return Err(crate::error::value("pairs must have shape (n, 2)"));
    }
    // `as_slice` also accepts Fortran order, which would silently transpose the rows.
    let flat = pairs
        .as_slice()
        .ok()
        .filter(|_| pairs.is_c_contiguous())
        .ok_or_else(|| {
            crate::error::value("pairs must be C-contiguous; call numpy.ascontiguousarray")
        })?;
    flat.as_chunks::<2>()
        .0
        .iter()
        .map(|&[a, b]| {
            let index = |value: u64| {
                usize::try_from(value).map_err(|_| crate::error::value("pair index exceeds usize"))
            };
            Ok((index(a)?, index(b)?))
        })
        .collect()
}

/// TM-scores of structure pairs, parallel across pairs.
///
/// With `pairs=None` every unordered pair is scored in condensed order
/// `(0,1), (0,2), ..., (n-2,n-1)`, structure `i` as the model and `j` as the
/// reference. Otherwise `pairs` is an `(m, 2)` unsigned array of structure
/// indices, `(model, reference)` per row.
#[pyfunction]
#[pyo3(signature = (structures, *, pairs=None, context=None))]
fn tm_scores<'py>(
    py: Python<'py>,
    structures: Vec<Readonly2<'py, f32>>,
    pairs: Option<Readonly2<'py, u64>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let views = views(&structures)?;
    let pairs = pairs.map(|array| pair_list(&array)).transpose()?;
    let scores = crate::execution::run(py, context, |context| match &pairs {
        Some(pairs) => molframe::compare::tm_score_pairs(&views, pairs, context),
        None => molframe::compare::tm_score_all_pairs(&views, context),
    })?
    .map_err(crate::error::kernel)?;
    // Release the NumPy read borrows now that the kernel is done.
    drop(views);
    drop(structures);
    Ok(scores.into_pyarray(py))
}

/// [`tm_scores`] behind an opt-in prefilter.
///
/// The exact stage drops a pair only when a proven upper bound on its
/// TM-score is below `threshold`, so every pair at or above it is scored
/// exactly as `tm_score` would. With `coarse_ratio` set, a heuristic stage also
/// drops pairs whose single full-length fit scores under
/// `coarse_ratio * threshold`; its recall is measured, not guaranteed. Dropped
/// pairs score `nan`.
#[pyfunction]
#[pyo3(signature = (structures, *, threshold, coarse_ratio=None, pairs=None, context=None))]
fn tm_screen<'py>(
    py: Python<'py>,
    structures: Vec<Readonly2<'py, f32>>,
    threshold: f64,
    coarse_ratio: Option<f64>,
    pairs: Option<Readonly2<'py, u64>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyTmScreen> {
    let views = views(&structures)?;
    let pairs = pairs.map(|array| pair_list(&array)).transpose()?;
    let mut filter =
        molframe::compare::TmPrefilter::new(threshold).map_err(crate::error::kernel)?;
    if let Some(ratio) = coarse_ratio {
        filter = filter
            .with_coarse_ratio(ratio)
            .map_err(crate::error::kernel)?;
    }
    let screen = crate::execution::run(py, context, |context| match &pairs {
        Some(pairs) => molframe::compare::tm_screen_pairs(&views, pairs, &filter, context),
        None => molframe::compare::tm_screen_all_pairs(&views, &filter, context),
    })?
    .map_err(crate::error::kernel)?;
    drop(views);
    drop(structures);
    Ok(PyTmScreen {
        scores: screen.scores.into_pyarray(py).unbind(),
        evaluated: screen.evaluated,
        pruned_bound: screen.pruned_bound,
        pruned_coarse: screen.pruned_coarse,
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTmScreen>()?;
    module.add_function(wrap_pyfunction!(tm_scores, module)?)?;
    module.add_function(wrap_pyfunction!(tm_screen, module)?)
}
