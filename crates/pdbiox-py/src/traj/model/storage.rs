//! Owned and retained-`NumPy` trajectory storage with explicit materialization.

use super::super::arrays;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyReadonlyArray3, PyUntypedArrayMethods,
};
use pdbiox::UnitCell;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::ops::Deref;
use std::sync::Arc;

const VECTOR_WIDTH: usize = 3;
const CELL_WIDTH: usize = 6;
const CONTIGUOUS: &str =
    "trajectory arrays must be C-contiguous; call numpy.ascontiguousarray explicitly";
const STALE_SHAPE: &str = "trajectory array shape changed after the view was created";

#[derive(Clone)]
pub(super) enum VectorStorage {
    Owned(Arc<[f32]>),
    Retained(Arc<Py<PyArray3<f32>>>),
}

#[derive(Clone)]
pub(super) enum ScalarStorage {
    Owned(Arc<[f64]>),
    Retained(Arc<Py<PyArray1<f64>>>),
}

#[derive(Clone)]
pub(super) enum CellStorage {
    Owned(Arc<[f64]>),
    Retained(Arc<Py<PyArray2<f64>>>),
}

#[derive(Clone)]
pub(super) enum IntegerStorage {
    Owned(Arc<[i64]>),
    Retained(Arc<Py<PyArray1<i64>>>),
}

impl VectorStorage {
    pub(super) fn from_array(values: PyReadonlyArray3<'_, f32>, copy: bool) -> PyResult<Self> {
        let contiguous_values = contiguous(&values)?;
        if copy {
            Ok(Self::Owned(Arc::from(contiguous_values.to_vec())))
        } else {
            Ok(Self::Retained(Arc::new(values.deref().clone().unbind())))
        }
    }

    pub(super) const fn owned(values: Arc<[f32]>) -> Self {
        Self::Owned(values)
    }

    pub(super) fn is_retained(&self) -> bool {
        matches!(self, Self::Retained(_))
    }

    pub(super) fn array<'py>(
        &self,
        py: Python<'py>,
        frames: usize,
        atoms: usize,
    ) -> PyResult<Bound<'py, PyArray3<f32>>> {
        match self {
            Self::Owned(values) => arrays::vectors(py, Arc::clone(values), frames, atoms),
            Self::Retained(values) => {
                let array = values.bind(py);
                vector_shape(array, frames, atoms)?;
                contiguous_array(array)?;
                Ok(array.clone())
            }
        }
    }

    pub(super) fn frame(
        &self,
        py: Python<'_>,
        frame: usize,
        frames: usize,
        atoms: usize,
    ) -> PyResult<Vec<[f32; VECTOR_WIDTH]>> {
        let values = match self {
            Self::Owned(values) => values.as_ref(),
            Self::Retained(values) => {
                let array = values.bind(py);
                vector_shape(array, frames, atoms)?;
                contiguous_array(array)?;
                let borrowed = array.readonly();
                return vector_frame(
                    borrowed
                        .as_slice()
                        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?,
                    frame,
                    atoms,
                );
            }
        };
        vector_frame(values, frame, atoms)
    }
}

impl ScalarStorage {
    pub(super) fn from_array(values: PyReadonlyArray1<'_, f64>, copy: bool) -> PyResult<Self> {
        let contiguous_values = contiguous(&values)?;
        if copy {
            Ok(Self::Owned(Arc::from(contiguous_values.to_vec())))
        } else {
            Ok(Self::Retained(Arc::new(values.deref().clone().unbind())))
        }
    }

    pub(super) const fn owned(values: Arc<[f64]>) -> Self {
        Self::Owned(values)
    }

    pub(super) fn array<'py>(
        &self,
        py: Python<'py>,
        frames: usize,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        match self {
            Self::Owned(values) => arrays::scalars(py, Arc::clone(values)),
            Self::Retained(values) => {
                let array = values.bind(py);
                scalar_shape(array, frames)?;
                contiguous_array(array)?;
                Ok(array.clone())
            }
        }
    }

    pub(super) fn value(&self, py: Python<'_>, index: usize, frames: usize) -> PyResult<f64> {
        match self {
            Self::Owned(values) => values
                .get(index)
                .copied()
                .ok_or_else(|| PyValueError::new_err(STALE_SHAPE)),
            Self::Retained(values) => {
                let array = values.bind(py);
                scalar_shape(array, frames)?;
                contiguous_array(array)?;
                array
                    .readonly()
                    .get(index)
                    .copied()
                    .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))
            }
        }
    }
}

impl CellStorage {
    pub(super) fn from_array(values: PyReadonlyArray2<'_, f64>, copy: bool) -> PyResult<Self> {
        let contiguous_values = contiguous(&values)?;
        if copy {
            Ok(Self::Owned(Arc::from(contiguous_values.to_vec())))
        } else {
            Ok(Self::Retained(Arc::new(values.deref().clone().unbind())))
        }
    }

    pub(super) const fn owned(values: Arc<[f64]>) -> Self {
        Self::Owned(values)
    }

    pub(super) fn array<'py>(
        &self,
        py: Python<'py>,
        frames: usize,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        match self {
            Self::Owned(values) => arrays::cells(py, Arc::clone(values), frames),
            Self::Retained(values) => {
                let array = values.bind(py);
                cell_shape(array, frames)?;
                contiguous_array(array)?;
                Ok(array.clone())
            }
        }
    }

    pub(super) fn cell(&self, py: Python<'_>, index: usize, frames: usize) -> PyResult<UnitCell> {
        let values = match self {
            Self::Owned(values) => values.as_ref(),
            Self::Retained(values) => {
                let array = values.bind(py);
                cell_shape(array, frames)?;
                contiguous_array(array)?;
                let borrowed = array.readonly();
                return cell_from_values(
                    borrowed
                        .as_slice()
                        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?,
                    index,
                );
            }
        };
        cell_from_values(values, index)
    }
}

impl IntegerStorage {
    pub(super) fn from_array(values: PyReadonlyArray1<'_, i64>, copy: bool) -> PyResult<Self> {
        let contiguous_values = contiguous(&values)?;
        if copy {
            Ok(Self::Owned(Arc::from(contiguous_values.to_vec())))
        } else {
            Ok(Self::Retained(Arc::new(values.deref().clone().unbind())))
        }
    }

    pub(super) const fn owned(values: Arc<[i64]>) -> Self {
        Self::Owned(values)
    }

    pub(super) fn array<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<i64>>> {
        match self {
            Self::Owned(values) => arrays::integers(py, Arc::clone(values)),
            Self::Retained(values) => {
                let array = values.bind(py);
                contiguous_array(array)?;
                Ok(array.clone())
            }
        }
    }

    pub(super) fn values(&self, py: Python<'_>, frames: usize) -> PyResult<Vec<i64>> {
        match self {
            Self::Owned(values) if values.len() == frames => Ok(values.to_vec()),
            Self::Owned(_) => Err(PyValueError::new_err(STALE_SHAPE)),
            Self::Retained(values) => {
                let array = values.bind(py);
                scalar_shape(array, frames)?;
                contiguous_array(array)?;
                Ok(array
                    .readonly()
                    .as_slice()
                    .map_err(|_| PyValueError::new_err(CONTIGUOUS))?
                    .to_vec())
            }
        }
    }
}

fn contiguous<'py, T: numpy::Element, D: numpy::ndarray::Dimension>(
    values: &'py numpy::PyReadonlyArray<'py, T, D>,
) -> PyResult<&'py [T]> {
    values
        .as_slice()
        .map_err(|_| PyValueError::new_err(CONTIGUOUS))
}

fn contiguous_array<T: numpy::Element, D: numpy::ndarray::Dimension>(
    values: &Bound<'_, numpy::PyArray<T, D>>,
) -> PyResult<()> {
    values
        .is_c_contiguous()
        .then_some(())
        .ok_or_else(|| PyValueError::new_err(CONTIGUOUS))
}

fn vector_shape(values: &Bound<'_, PyArray3<f32>>, frames: usize, atoms: usize) -> PyResult<()> {
    (values.shape() == [frames, atoms, VECTOR_WIDTH])
        .then_some(())
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))
}

fn scalar_shape<T: numpy::Element>(values: &Bound<'_, PyArray1<T>>, frames: usize) -> PyResult<()> {
    (values.shape() == [frames])
        .then_some(())
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))
}

fn cell_shape(values: &Bound<'_, PyArray2<f64>>, frames: usize) -> PyResult<()> {
    (values.shape() == [frames, CELL_WIDTH])
        .then_some(())
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))
}

fn vector_frame(values: &[f32], frame: usize, atoms: usize) -> PyResult<Vec<[f32; VECTOR_WIDTH]>> {
    let width = atoms
        .checked_mul(VECTOR_WIDTH)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    let start = frame
        .checked_mul(width)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    let end = start
        .checked_add(width)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    let frame = values
        .get(start..end)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    Ok(frame
        .chunks_exact(VECTOR_WIDTH)
        .map(|value| [value[0], value[1], value[2]])
        .collect())
}

fn cell_from_values(values: &[f64], frame: usize) -> PyResult<UnitCell> {
    let start = frame
        .checked_mul(CELL_WIDTH)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    let end = start
        .checked_add(CELL_WIDTH)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    let cell = values
        .get(start..end)
        .ok_or_else(|| PyValueError::new_err(STALE_SHAPE))?;
    Ok(UnitCell {
        lengths: [cell[0], cell[1], cell[2]],
        angles: [cell[3], cell[4], cell[5]],
    })
}
