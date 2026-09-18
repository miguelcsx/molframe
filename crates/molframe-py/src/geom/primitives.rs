//! Vectorized bindings for primitive geometric kernels and fluctuations.

use super::arrays::{distance_matrix_value, matrix_error};
use super::borrowed_coordinates;
use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyReadonlyArray3, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn distance_matrix_between<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f32>,
    right: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let left = borrowed_coordinates(&left)?;
    let right = borrowed_coordinates(&right)?;
    let matrix = py
        .detach(move || molframe::distance_matrix_between(left, right))
        .map_err(matrix_error)?;
    distance_matrix_value(py, matrix)
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
    let values = frames.as_slice().map_err(|_| {
        PyValueError::new_err(
            "frames must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let (values, remainder) = values.as_chunks::<3>();
    if !remainder.is_empty() {
        return Err(PyValueError::new_err(
            "frames must have shape (frames, atoms, 3)",
        ));
    }
    let borrowed = if atom_count == 0 {
        std::iter::repeat_n(&values[..0], shape[0]).collect::<Vec<_>>()
    } else {
        values.chunks_exact(atom_count).collect::<Vec<_>>()
    };
    let values = py
        .detach(|| molframe::rmsf(&borrowed))
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
            let left = borrowed_coordinates(&left)?.to_vec();
            let right = borrowed_coordinates(&right)?.to_vec();
            equal_rows(&left, &right)?;
            let values = py.detach(|| {
                left.iter()
                    .zip(right.iter())
                    .map(|(&left, &right)| $kernel(left, right))
                    .collect::<Vec<_>>()
            });
            Ok(Array1::from_vec(values).into_pyarray(py))
        }
    };
}

binary_scalar!(distance_squared, molframe::distance_squared);

#[pyfunction]
pub(crate) fn distance<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f32>,
    right: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let left = borrowed_coordinates(&left)?.to_vec();
    let right = borrowed_coordinates(&right)?.to_vec();
    equal_rows(&left, &right)?;
    let values = py.detach(move || {
        let mut output = vec![0.0; left.len()];
        molframe::distances_into(&left, &right, &mut output)
            .map(|_| output)
            .map_err(batch_error)
    })?;
    Ok(readonly_array(py, values))
}

#[pyfunction]
#[pyo3(signature = (left, right))]
pub(crate) fn displacement<'py>(
    py: Python<'py>,
    left: PyReadonlyArray2<'_, f32>,
    right: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let left = borrowed_coordinates(&left)?;
    let right = borrowed_coordinates(&right)?;
    equal_rows(left, right)?;
    let values = py.detach(|| {
        left.iter()
            .zip(right.iter())
            .map(|(&from, &to)| molframe::displacement(from, to))
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
    let left = borrowed_vectors64(&left)?;
    let right = borrowed_vectors64(&right)?;
    equal_rows(left, right)?;
    let values = py.detach(|| {
        left.iter()
            .zip(right.iter())
            .map(|(&left, &right)| molframe::dot(left, right))
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
    let left = borrowed_vectors64(&left)?;
    let right = borrowed_vectors64(&right)?;
    equal_rows(left, right)?;
    let values = py.detach(|| {
        left.iter()
            .zip(right.iter())
            .map(|(&left, &right)| molframe::cross(left, right))
            .collect::<Vec<_>>()
    });
    triples(py, values)
}

#[pyfunction]
pub(crate) fn norm<'py>(
    py: Python<'py>,
    vectors: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let vectors = borrowed_vectors64(&vectors)?;
    let values = py.detach(|| {
        vectors
            .iter()
            .map(|&vector| molframe::norm(vector))
            .collect::<Vec<_>>()
    });
    Ok(Array1::from_vec(values).into_pyarray(py))
}

#[pyfunction]
pub(crate) fn normalise(
    py: Python<'_>,
    vectors: PyReadonlyArray2<'_, f64>,
) -> PyResult<Vec<Option<[f64; 3]>>> {
    let vectors = borrowed_vectors64(&vectors)?;
    Ok(py.detach(|| {
        vectors
            .iter()
            .map(|&vector| molframe::normalise(vector))
            .collect()
    }))
}

#[pyfunction]
pub(crate) fn angle<'py>(
    py: Python<'py>,
    first: PyReadonlyArray2<'_, f32>,
    vertex: PyReadonlyArray2<'_, f32>,
    third: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let first = borrowed_coordinates(&first)?.to_vec();
    let vertex = borrowed_coordinates(&vertex)?.to_vec();
    let third = borrowed_coordinates(&third)?.to_vec();
    equal_rows(&first, &vertex)?;
    equal_rows(&first, &third)?;
    let values = py.detach(move || {
        let mut output = vec![0.0; first.len()];
        molframe::angles_into(&first, &vertex, &third, &mut output)
            .map(|_| output)
            .map_err(batch_error)
    })?;
    Ok(readonly_array(py, values))
}

#[pyfunction]
pub(crate) fn dihedral<'py>(
    py: Python<'py>,
    first: PyReadonlyArray2<'_, f32>,
    second: PyReadonlyArray2<'_, f32>,
    third: PyReadonlyArray2<'_, f32>,
    fourth: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let first = borrowed_coordinates(&first)?.to_vec();
    let second = borrowed_coordinates(&second)?.to_vec();
    let third = borrowed_coordinates(&third)?.to_vec();
    let fourth = borrowed_coordinates(&fourth)?.to_vec();
    equal_rows(&first, &second)?;
    equal_rows(&first, &third)?;
    equal_rows(&first, &fourth)?;
    let values = py.detach(move || {
        let mut output = vec![0.0; first.len()];
        molframe::torsions_into(&first, &second, &third, &fourth, &mut output)
            .map(|_| output)
            .map_err(batch_error)
    })?;
    Ok(readonly_array(py, values))
}

#[pyfunction]
pub(crate) fn degrees<'py>(
    py: Python<'py>,
    radians: PyReadonlyArray1<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let values = radians.as_slice()?;
    let values = py.detach(|| {
        values
            .iter()
            .map(|&value| molframe::degrees(value))
            .collect::<Vec<_>>()
    });
    Ok(Array1::from_vec(values).into_pyarray(py))
}

fn borrowed_vectors64<'a>(values: &'a PyReadonlyArray2<'_, f64>) -> PyResult<&'a [[f64; 3]]> {
    if values.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("vectors must have shape (n, 3)"));
    }
    let values = values.as_slice().map_err(|_| {
        PyValueError::new_err(
            "vectors must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let (vectors, remainder) = values.as_chunks::<3>();
    if remainder.is_empty() {
        Ok(vectors)
    } else {
        Err(PyValueError::new_err("vectors must have shape (n, 3)"))
    }
}

fn equal_rows<T>(left: &[T], right: &[T]) -> PyResult<()> {
    if left.len() == right.len() {
        Ok(())
    } else {
        Err(PyValueError::new_err("vector batches differ in length"))
    }
}

fn readonly_array<T: numpy::Element>(py: Python<'_>, values: Vec<T>) -> Bound<'_, PyArray1<T>> {
    let array = Array1::from_vec(values).into_pyarray(py);
    array.readwrite().make_nonwriteable();
    array
}

fn batch_error(error: molframe::BatchGeometryError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn triples(py: Python<'_>, values: Vec<[f64; 3]>) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let rows = values.len();
    Array2::from_shape_vec((rows, 3), values.into_iter().flatten().collect())
        .map(|array| array.into_pyarray(py))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
