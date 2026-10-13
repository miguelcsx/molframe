//! Stateless geometry over `NumPy` coordinate arrays.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;

pub(crate) fn coordinates<'a>(array: &'a PyReadonlyArray2<'_, f32>) -> PyResult<&'a [[f32; 3]]> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(crate::error::value("coordinates must have shape (n, 3)"));
    }
    let contiguous = array.as_slice().map_err(|_| {
        crate::error::value("coordinates must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    Ok(contiguous.as_chunks::<3>().0)
}

#[pyfunction]
pub(crate) fn centroid(array: &Bound<'_, PyArray2<f32>>) -> PyResult<Option<[f64; 3]>> {
    let array = array.readonly();
    Ok(molframe::geometry::centroid(coordinates(&array)?))
}

#[pyfunction]
pub(crate) fn rmsd(
    mobile: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    let mobile = mobile.readonly();
    let reference = reference.readonly();
    molframe::geometry::rmsd(coordinates(&mobile)?, coordinates(&reference)?)
        .map_err(crate::error::kernel)
}

#[derive(Clone, Copy)]
enum Width {
    Double,
    Single,
}

fn width(dtype: &str) -> PyResult<Width> {
    match dtype {
        "float64" | "f64" => Ok(Width::Double),
        "float32" | "f32" => Ok(Width::Single),
        other => Err(crate::error::value(format!(
            "dtype must be \"float64\" or \"float32\", not {other:?}"
        ))),
    }
}

/// All pairwise distances, charged to the context memory budget.
///
/// `dtype="float32"` halves the footprint; `condensed=True` returns the upper
/// triangle as a 1-D array in `scipy.spatial.distance.pdist` order.
#[pyfunction]
#[pyo3(signature = (array, *, dtype="float64", condensed=false, context=None))]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    array: &Bound<'py, PyArray2<f32>>,
    dtype: &str,
    condensed: bool,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyAny>> {
    use molframe::geometry as geom;
    let width = width(dtype)?;
    let array = array.readonly();
    let positions = coordinates(&array)?;
    let size = positions.len();
    let run =
        |work: &(
              dyn Fn(&molframe::ExecutionContext) -> Result<Output, geom::MatrixError> + Sync
          )| { crate::execution::run(py, context, work)?.map_err(crate::error::kernel) };
    let output = match (width, condensed) {
        (Width::Double, false) => run(&|context| {
            geom::distance_matrix_with_context(positions, context)
                .map(|matrix| Output::Double(matrix.into_values()))
        })?,
        (Width::Single, false) => run(&|context| {
            geom::distance_matrix_f32_with_context(positions, context).map(Output::Single)
        })?,
        (Width::Double, true) => run(&|context| {
            geom::condensed_distances_with_context(positions, context).map(Output::Double)
        })?,
        (Width::Single, true) => run(&|context| {
            geom::condensed_distances_f32_with_context(positions, context).map(Output::Single)
        })?,
    };
    output.into_python(py, (!condensed).then_some((size, size)))
}

/// One block of `distance_matrix` rows `[start, stop)`, as `(stop - start, n)`.
#[pyfunction]
#[pyo3(name = "_distance_rows", signature = (array, start, stop, *, dtype="float64", context=None))]
pub(crate) fn distance_rows<'py>(
    py: Python<'py>,
    array: &Bound<'py, PyArray2<f32>>,
    start: usize,
    stop: usize,
    dtype: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyAny>> {
    use molframe::geometry as geom;
    let width = width(dtype)?;
    let array = array.readonly();
    let positions = coordinates(&array)?;
    let size = positions.len();
    if start > stop || stop > size {
        return Err(crate::error::value(format!(
            "row range {start}..{stop} is outside a matrix of {size} rows"
        )));
    }
    let rows = stop - start;
    let element = match width {
        Width::Double => size_of::<f64>(),
        Width::Single => size_of::<f32>(),
    };
    let output = crate::execution::run(
        py,
        context,
        |context| -> Result<Output, geom::MatrixError> {
            let bytes = rows.saturating_mul(size).saturating_mul(element);
            let _held =
                context
                    .try_reserve(bytes)
                    .map_err(|_| geom::MatrixError::MemoryBudget {
                        requested: bytes,
                        available: context.available_bytes(),
                    })?;
            match width {
                Width::Double => {
                    let mut block = vec![0.0_f64; rows * size];
                    geom::distance_rows_into(positions, start..stop, &mut block)?;
                    Ok(Output::Double(block))
                }
                Width::Single => {
                    let mut block = vec![0.0_f32; rows * size];
                    geom::distance_rows_f32_into(positions, start..stop, &mut block)?;
                    Ok(Output::Single(block))
                }
            }
        },
    )?
    .map_err(crate::error::kernel)?;
    output.into_python(py, Some((rows, size)))
}

/// Rows per block that keep a stream of `distance_rows` blocks within budget.
#[pyfunction]
#[pyo3(name = "_distance_block_rows", signature = (array, *, dtype="float64", context=None))]
pub(crate) fn distance_block_rows(
    py: Python<'_>,
    array: &Bound<'_, PyArray2<f32>>,
    dtype: &str,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<usize> {
    let width = width(dtype)?;
    let size = coordinates(&array.readonly())?.len();
    let element = match width {
        Width::Double => size_of::<f64>(),
        Width::Single => size_of::<f32>(),
    };
    crate::execution::run(py, context, |context| {
        molframe::geometry::distance_block_rows(size, element, context)
    })?
    .map_err(crate::error::kernel)
}

enum Output {
    Double(Vec<f64>),
    Single(Vec<f32>),
}

impl Output {
    fn into_python(
        self,
        py: Python<'_>,
        shape: Option<(usize, usize)>,
    ) -> PyResult<Bound<'_, PyAny>> {
        match (self, shape) {
            (Self::Double(values), Some(shape)) => {
                Ok(values.into_pyarray(py).reshape(shape)?.into_any())
            }
            (Self::Single(values), Some(shape)) => {
                Ok(values.into_pyarray(py).reshape(shape)?.into_any())
            }
            (Self::Double(values), None) => Ok(values.into_pyarray(py).into_any()),
            (Self::Single(values), None) => Ok(values.into_pyarray(py).into_any()),
        }
    }
}
