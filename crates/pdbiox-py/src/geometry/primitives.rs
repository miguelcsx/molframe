//! Vectorized bindings for primitive geometric kernels and fluctuations.

use super::arrays::matrix_error;
use super::coordinates;
use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn distance_matrix_between<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f32>,
    right: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let left = coordinates(left)?;
    let right = coordinates(right)?;
    let matrix = py
        .detach(move || pdbiox::distance_matrix_between(&left, &right))
        .map_err(matrix_error)?;
    matrix_array(py, &matrix)
}

#[pyfunction]
pub(crate) fn rmsf<'py>(
    py: Python<'py>,
    frames: PyReadonlyArray3<'_, f32>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let shape = frames.shape();
    if shape.len() != 3 || shape[2] != 3 {
        return Err(PyValueError::new_err(
            "frames must have shape (frames, atoms, 3)",
        ));
    }
    let atom_count = shape[1];
    let owned = frames
        .as_array()
        .outer_iter()
        .map(|frame| {
            frame
                .outer_iter()
                .map(|row| [row[0], row[1], row[2]])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    drop(frames);
    let values = py
        .detach(|| {
            let borrowed = owned.iter().map(Vec::as_slice).collect::<Vec<_>>();
            pdbiox::rmsf(&borrowed)
        })
        .map_err(|error| PyValueError::new_err(format!("RMSF failed: {error:?}")))?;
    if values.len() != atom_count {
        return Err(PyValueError::new_err("RMSF returned an invalid shape"));
    }
    Ok(Array1::from_vec(values).into_pyarray(py))
}

macro_rules! binary_scalar {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name<'py>(
            py: Python<'py>,
            left: PyReadonlyArray2<'_, f32>,
            right: PyReadonlyArray2<'_, f32>,
        ) -> PyResult<Bound<'py, PyArray1<f64>>> {
            let left = coordinates(left)?;
            let right = coordinates(right)?;
            equal_rows(&left, &right)?;
            let values = py.detach(|| {
                left.into_iter()
                    .zip(right)
                    .map(|(left, right)| $kernel(left, right))
                    .collect::<Vec<_>>()
            });
            Ok(Array1::from_vec(values).into_pyarray(py))
        }
    };
}

binary_scalar!(distance, pdbiox::distance);
binary_scalar!(distance_squared, pdbiox::distance_squared);

#[pyfunction]
pub(crate) fn displacement<'py>(
    py: Python<'py>,
    from: PyReadonlyArray2<'_, f32>,
    to: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let from = coordinates(from)?;
    let to = coordinates(to)?;
    equal_rows(&from, &to)?;
    let values = py.detach(|| {
        from.into_iter()
            .zip(to)
            .map(|(from, to)| pdbiox::displacement(from, to))
            .collect::<Vec<_>>()
    });
    triples(py, values)
}

#[pyfunction]
pub(crate) fn dot<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f64>,
    right: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let left = vectors64(left)?;
    let right = vectors64(right)?;
    equal_rows(&left, &right)?;
    let values = py.detach(|| {
        left.into_iter()
            .zip(right)
            .map(|(left, right)| pdbiox::dot(left, right))
            .collect::<Vec<_>>()
    });
    Ok(Array1::from_vec(values).into_pyarray(py))
}

#[pyfunction]
pub(crate) fn cross<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f64>,
    right: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let left = vectors64(left)?;
    let right = vectors64(right)?;
    equal_rows(&left, &right)?;
    let values = py.detach(|| {
        left.into_iter()
            .zip(right)
            .map(|(left, right)| pdbiox::cross(left, right))
            .collect::<Vec<_>>()
    });
    triples(py, values)
}

#[pyfunction]
pub(crate) fn norm<'py>(
    py: Python<'py>,
    vectors: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let vectors = vectors64(vectors)?;
    let values = py.detach(|| vectors.into_iter().map(pdbiox::norm).collect::<Vec<_>>());
    Ok(Array1::from_vec(values).into_pyarray(py))
}

#[pyfunction]
pub(crate) fn normalise(
    py: Python<'_>,
    vectors: PyReadonlyArray2<'_, f64>,
) -> PyResult<Vec<Option<[f64; 3]>>> {
    let vectors = vectors64(vectors)?;
    Ok(py.detach(|| vectors.into_iter().map(pdbiox::normalise).collect()))
}

#[pyfunction]
pub(crate) fn angle(
    py: Python<'_>,
    first: PyReadonlyArray2<'_, f32>,
    vertex: PyReadonlyArray2<'_, f32>,
    third: PyReadonlyArray2<'_, f32>,
) -> PyResult<Vec<Option<f64>>> {
    let first = coordinates(first)?;
    let vertex = coordinates(vertex)?;
    let third = coordinates(third)?;
    equal_rows(&first, &vertex)?;
    equal_rows(&first, &third)?;
    Ok(py.detach(|| {
        first
            .into_iter()
            .zip(vertex)
            .zip(third)
            .map(|((first, vertex), third)| pdbiox::angle(first, vertex, third))
            .collect()
    }))
}

#[pyfunction]
pub(crate) fn dihedral(
    py: Python<'_>,
    first: PyReadonlyArray2<'_, f32>,
    second: PyReadonlyArray2<'_, f32>,
    third: PyReadonlyArray2<'_, f32>,
    fourth: PyReadonlyArray2<'_, f32>,
) -> PyResult<Vec<Option<f64>>> {
    let first = coordinates(first)?;
    let second = coordinates(second)?;
    let third = coordinates(third)?;
    let fourth = coordinates(fourth)?;
    equal_rows(&first, &second)?;
    equal_rows(&first, &third)?;
    equal_rows(&first, &fourth)?;
    Ok(py.detach(|| {
        first
            .into_iter()
            .zip(second)
            .zip(third)
            .zip(fourth)
            .map(|(((a, b), c), d)| pdbiox::dihedral(a, b, c, d))
            .collect()
    }))
}

#[pyfunction]
pub(crate) fn degrees<'py>(
    py: Python<'py>,
    radians: PyReadonlyArray1<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let values = radians.as_slice()?.to_vec();
    drop(radians);
    let values = py.detach(|| values.into_iter().map(pdbiox::degrees).collect::<Vec<_>>());
    Ok(Array1::from_vec(values).into_pyarray(py))
}

fn vectors64(values: PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    if values.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("vectors must have shape (n, 3)"));
    }
    let output = values
        .as_array()
        .outer_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    drop(values);
    Ok(output)
}

fn equal_rows<T>(left: &[T], right: &[T]) -> PyResult<()> {
    if left.len() == right.len() {
        Ok(())
    } else {
        Err(PyValueError::new_err("vector batches differ in length"))
    }
}

fn triples(py: Python<'_>, values: Vec<[f64; 3]>) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let rows = values.len();
    Array2::from_shape_vec((rows, 3), values.into_iter().flatten().collect())
        .map(|array| array.into_pyarray(py))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

fn matrix_array<'py>(
    py: Python<'py>,
    matrix: &pdbiox::DistanceMatrix,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    Array2::from_shape_vec(
        (matrix.rows(), matrix.columns()),
        matrix.as_slice().to_vec(),
    )
    .map(|array| array.into_pyarray(py))
    .map_err(|error| PyValueError::new_err(error.to_string()))
}
