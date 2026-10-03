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
#[pyo3(signature = (coordinates, cutoff, *, backend="auto"))]
fn neighbor_pairs<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    cutoff: f32,
    backend: &str,
) -> PyResult<Pairs<'py>> {
    let backend = crate::backend::parse(backend)?;
    let array = coordinates.readonly();
    let positions = crate::bindings::coordinates(&array)?;
    let everything = molframe::engine::core::AtomSelection::All(
        u32::try_from(positions.len())
            .map_err(|_| crate::error::value("too many atoms for a 32-bit index"))?,
    );
    let pairs = py
        .detach(|| {
            molframe::spatial::pairs_within(
                positions,
                &everything,
                &everything,
                cutoff,
                backend,
                None,
                &molframe::ExecutionContext::default(),
            )
        })
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
    module.add_function(wrap_pyfunction!(neighbor_pairs, module)?)
}
