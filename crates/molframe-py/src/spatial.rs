//! Mechanical adapters for fixed-radius neighbour search over coordinates.

use numpy::{PyArray1, PyArray2, PyArrayMethods, ToPyArray};
use pyo3::prelude::*;

type Pairs<'py> = (
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<f32>>,
);

/// Unique unordered pairs no further apart than `cutoff`, sorted by index.
///
/// Returns `(first, second, distance)` arrays with `first < second`.
#[pyfunction]
#[pyo3(signature = (coordinates, cutoff, *, backend="auto", context=None))]
fn neighbor_pairs<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    cutoff: f32,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Pairs<'py>> {
    let backend = crate::backend::parse(backend)?;
    let array = coordinates.readonly();
    let positions = crate::bindings::coordinates(&array)?;
    let everything = molframe::engine::core::AtomSelection::All(
        u32::try_from(positions.len())
            .map_err(|_| crate::error::value("too many atoms for a 32-bit index"))?,
    );
    let pairs = crate::execution::run(py, context, |context| {
        molframe::spatial::pairs_within(
            positions,
            &everything,
            &everything,
            cutoff,
            backend,
            None,
            context,
        )
    })?
    .map_err(crate::error::kernel)?;
    let first: Vec<u32> = pairs.iter().map(|pair| pair.first).collect();
    let second: Vec<u32> = pairs.iter().map(|pair| pair.second).collect();
    let distance: Vec<f32> = pairs
        .iter()
        .map(|pair| pair.distance_squared.sqrt())
        .collect();
    Ok((
        first.to_pyarray(py),
        second.to_pyarray(py),
        distance.to_pyarray(py),
    ))
}

/// Cross-set pairs only, avoiding the quadratic within-set output of dense clouds.
#[pyfunction]
#[pyo3(signature = (first, second, cutoff, *, backend="auto", context=None))]
fn cross_pairs<'py>(
    py: Python<'py>,
    first: &Bound<'py, PyArray2<f32>>,
    second: &Bound<'py, PyArray2<f32>>,
    cutoff: f32,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Pairs<'py>> {
    let backend = crate::backend::parse(backend)?;
    let first = first.readonly();
    let second = second.readonly();
    let first = crate::bindings::coordinates(&first)?;
    let second = crate::bindings::coordinates(&second)?;
    let pairs = crate::execution::run(py, context, |context| {
        molframe::spatial::cross_pairs(first, second, cutoff, backend, context)
    })?
    .map_err(crate::error::kernel)?;
    let first: Vec<u32> = pairs.iter().map(|pair| pair.first).collect();
    let second: Vec<u32> = pairs.iter().map(|pair| pair.second).collect();
    let distance: Vec<f32> = pairs
        .iter()
        .map(|pair| pair.distance_squared.sqrt())
        .collect();
    Ok((
        first.to_pyarray(py),
        second.to_pyarray(py),
        distance.to_pyarray(py),
    ))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(neighbor_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(cross_pairs, module)?)
}
