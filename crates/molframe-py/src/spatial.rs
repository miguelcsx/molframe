//! Mechanical adapters for fixed-radius neighbour search over coordinates.

use crate::bindings::array::Array2;
use crate::crystal::PyUnitCell;
use molframe::spatial::{KnnOptions, NeighborPair, PairColumns};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

type Pairs<'py> = (
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<f32>>,
);

fn periodic_box(cell: Option<&PyUnitCell>) -> PyResult<Option<molframe::spatial::PeriodicBox>> {
    cell.map(|cell| {
        molframe::spatial::PeriodicBox::from_cell(cell.parameters()).map_err(crate::error::kernel)
    })
    .transpose()
}

/// Struct-of-arrays columns handed to `NumPy` without a further copy.
fn columns(py: Python<'_>, columns: PairColumns) -> Pairs<'_> {
    (
        columns.first.into_pyarray(py),
        columns.second.into_pyarray(py),
        columns.distance.into_pyarray(py),
    )
}

/// Unique unordered pairs no further apart than `cutoff`, sorted by index.
///
/// Returns `(first, second, distance)` arrays with `first < second`. With a
/// `cell` distances are minimum-image; without one they are plain Euclidean.
#[pyfunction]
#[pyo3(signature = (coordinates, cutoff, *, cell=None, backend="auto", context=None))]
fn neighbor_pairs<'py>(
    py: Python<'py>,
    coordinates: Array2<'_, 'py, f32>,
    cutoff: f32,
    cell: Option<&PyUnitCell>,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Pairs<'py>> {
    let backend = crate::backend::parse(backend)?;
    let periodic = periodic_box(cell)?;
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
            periodic.as_ref(),
            context,
        )
    })?
    .map_err(crate::error::kernel)?;
    Ok(columns(py, PairColumns::from_pairs(&pairs)))
}

/// Cross-set pairs only, avoiding the quadratic within-set output of dense clouds.
/// With a `cell` distances are minimum-image.
#[pyfunction]
#[pyo3(signature = (first, second, cutoff, *, cell=None, backend="auto", context=None))]
fn cross_pairs<'py>(
    py: Python<'py>,
    first: Array2<'_, 'py, f32>,
    second: Array2<'_, 'py, f32>,
    cutoff: f32,
    cell: Option<&PyUnitCell>,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Pairs<'py>> {
    let backend = crate::backend::parse(backend)?;
    let periodic = periodic_box(cell)?;
    let first = first.readonly();
    let second = second.readonly();
    let first = crate::bindings::coordinates(&first)?;
    let second = crate::bindings::coordinates(&second)?;
    let pairs = crate::execution::run(py, context, |context| {
        molframe::spatial::cross_pairs_periodic(
            first,
            second,
            cutoff,
            backend,
            periodic.as_ref(),
            context,
        )
    })?
    .map_err(crate::error::kernel)?;
    Ok(columns(py, PairColumns::from_cross(&pairs)))
}

/// Iterator over `(first, second, distance)` chunks of a pair list.
#[pyclass(name = "PairChunks", module = "molframe.spatial", skip_from_py_object)]
struct PairChunks {
    pairs: Vec<NeighborPair>,
    cursor: usize,
    size: usize,
}

#[pymethods]
impl PairChunks {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> Option<Pairs<'py>> {
        if self.cursor >= self.pairs.len() {
            self.pairs = Vec::new();
            return None;
        }
        let end = self.pairs.len().min(self.cursor + self.size);
        let chunk = self.pairs.get(self.cursor..end)?;
        let out = columns(py, PairColumns::from_pairs(chunk));
        self.cursor = end;
        Some(out)
    }
}

/// Like `neighbor_pairs`, yielding bounded `(first, second, distance)` chunks.
///
/// The pair list is built once in 12 bytes per pair; each chunk then allocates
/// only `chunk_size` entries, so consumers that reduce chunks never hold the
/// three full-size `NumPy` columns as well.
#[pyfunction]
#[pyo3(signature = (coordinates, cutoff, *, chunk_size=1_048_576, cell=None, backend="auto", context=None))]
fn neighbor_pairs_chunks(
    py: Python<'_>,
    coordinates: Array2<'_, '_, f32>,
    cutoff: f32,
    chunk_size: usize,
    cell: Option<&PyUnitCell>,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PairChunks> {
    if chunk_size == 0 {
        return Err(crate::error::value("chunk_size must be positive"));
    }
    let backend = crate::backend::parse(backend)?;
    let periodic = periodic_box(cell)?;
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
            periodic.as_ref(),
            context,
        )
    })?
    .map_err(crate::error::kernel)?;
    Ok(PairChunks {
        pairs,
        cursor: 0,
        size: chunk_size,
    })
}

type Neighbors<'py> = (Bound<'py, PyArray2<u32>>, Bound<'py, PyArray2<f32>>);

/// The `k` nearest `points` of every query point, as `(indices, distances)`.
///
/// Both outputs have shape `(n_query, k)`, rows ordered by increasing distance
/// with ties by index. Missing neighbours (fewer than `k` points, or none within
/// `max_distance`) hold `4294967295` and `inf`. Without `query` the points query
/// themselves and `exclude_self` defaults to true; with `query` it defaults to
/// false. Distances are float32 and, with a `cell`, minimum-image.
#[pyfunction]
#[pyo3(signature = (points, k, *, query=None, max_distance=None, cell=None, exclude_self=None, backend="auto", context=None))]
fn knn<'py>(
    py: Python<'py>,
    points: Array2<'_, 'py, f32>,
    k: usize,
    query: Option<Array2<'_, 'py, f32>>,
    max_distance: Option<f32>,
    cell: Option<&PyUnitCell>,
    exclude_self: Option<bool>,
    backend: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Neighbors<'py>> {
    let backend = crate::backend::parse(backend)?;
    let periodic = periodic_box(cell)?;
    let array = points.readonly();
    let targets = crate::bindings::coordinates(&array)?;
    let options = KnnOptions::for_queries(max_distance, exclude_self, query.is_some());
    let result = if let Some(query) = query {
        let query = query.readonly();
        let query = crate::bindings::coordinates(&query)?;
        crate::execution::run(py, context, |context| {
            molframe::spatial::knn_between(
                targets,
                query,
                k,
                options,
                backend,
                periodic.as_ref(),
                context,
            )
        })?
    } else {
        let count = u32::try_from(targets.len())
            .map_err(|_| crate::error::value("too many atoms for a 32-bit index"))?;
        let all = molframe::engine::core::AtomSelection::All(count);
        crate::execution::run(py, context, |context| {
            molframe::spatial::knn_within(
                targets,
                &all,
                &all,
                k,
                options,
                backend,
                periodic.as_ref(),
                context,
            )
        })?
    }
    .map_err(crate::error::kernel)?;
    let rows = result.rows();
    let (indices, distances) = result.into_columns();
    let shape = [rows, k];
    Ok((
        indices.into_pyarray(py).reshape(shape)?,
        distances.into_pyarray(py).reshape(shape)?,
    ))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(neighbor_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(cross_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(neighbor_pairs_chunks, module)?)?;
    module.add_function(wrap_pyfunction!(knn, module)?)?;
    module.add_class::<PairChunks>()
}
