//! Vectorized density and surface kernels with caller-owned scientific policy.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::geometry::borrowed_coordinates;
use numpy::ndarray::ArrayView3;
use numpy::{PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "CartesianAxis", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCartesianAxis {
    X,
    Y,
    Z,
}

#[pyclass(name = "LinearDensityBin", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLinearDensityBin {
    #[pyo3(get)]
    lower: f32,
    #[pyo3(get)]
    upper: f32,
    #[pyo3(get)]
    weight: f64,
    #[pyo3(get)]
    density: f64,
}

#[pyclass(name = "LinearDensityOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyLinearDensityOptions {
    #[pyo3(get)]
    axis: PyCartesianAxis,
    #[pyo3(get)]
    minimum: f32,
    #[pyo3(get)]
    maximum: f32,
    #[pyo3(get)]
    bins: usize,
}

#[pymethods]
impl PyLinearDensityOptions {
    #[new]
    fn new(axis: PyCartesianAxis, minimum: f32, maximum: f32, bins: usize) -> Self {
        Self {
            axis,
            minimum,
            maximum,
            bins,
        }
    }
}

impl PyLinearDensityOptions {
    pub(crate) fn native(&self) -> pdbiox::analysis::LinearDensityOptions {
        pdbiox::analysis::LinearDensityOptions {
            axis: self.axis.into(),
            minimum: self.minimum,
            maximum: self.maximum,
            bins: self.bins,
        }
    }

    pub(crate) fn from_native(options: pdbiox::analysis::LinearDensityOptions) -> Self {
        Self {
            axis: options.axis.into(),
            minimum: options.minimum,
            maximum: options.maximum,
            bins: options.bins,
        }
    }
}

#[pyclass(name = "DensityGridSpec", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDensityGridSpec(pdbiox::analysis::DensityGridSpec);

#[pymethods]
impl PyDensityGridSpec {
    #[new]
    fn new(origin: [f32; 3], spacing: [f32; 3], shape: [usize; 3]) -> Self {
        Self(pdbiox::analysis::DensityGridSpec {
            origin,
            spacing,
            shape,
        })
    }

    #[getter]
    const fn origin(&self) -> [f32; 3] {
        self.0.origin
    }
    #[getter]
    const fn spacing(&self) -> [f32; 3] {
        self.0.spacing
    }
    #[getter]
    const fn shape(&self) -> [usize; 3] {
        self.0.shape
    }
}

impl PyDensityGridSpec {
    pub(crate) const fn native(&self) -> pdbiox::analysis::DensityGridSpec {
        self.0
    }

    pub(crate) const fn from_native(spec: pdbiox::analysis::DensityGridSpec) -> Self {
        Self(spec)
    }
}

#[pyclass(name = "DensityGrid", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDensityGrid {
    spec: pdbiox::analysis::DensityGridSpec,
    values: Vec<f64>,
    #[pyo3(get)]
    excluded_weight: f64,
}

#[pymethods]
impl PyDensityGrid {
    #[getter]
    const fn spec(&self) -> PyDensityGridSpec {
        PyDensityGridSpec(self.spec)
    }

    #[getter]
    fn density<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let (shape, length, pointer) = {
            let grid = slf.borrow();
            (grid.spec.shape, grid.values.len(), grid.values.as_ptr())
        };
        let expected = shape.iter().try_fold(1_usize, |product, dimension| {
            product.checked_mul(*dimension)
        });
        if expected != Some(length) {
            return Err(PyValueError::new_err(
                "density-grid buffer length does not match its declared shape",
            ));
        }
        // SAFETY: the checked grid shape covers the complete contiguous f64
        // allocation. `slf` is retained as the array base below, so the frozen
        // grid and its vector cannot be dropped while the view exists.
        let view = unsafe { ArrayView3::from_shape_ptr(shape, pointer) };
        // SAFETY: the frozen grid owns the immutable backing allocation for the
        // whole lifetime of the exported view, which is made read-only below.
        let density = unsafe { PyArray3::borrow_from_array(&view, slf.clone().into_any()) };
        let _readonly = density.readwrite().make_nonwriteable();
        Ok(density)
    }
}

#[pyfunction]
pub(crate) fn linear_density(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    weights: PyReadonlyArray1<'_, f64>,
    axis: PyCartesianAxis,
    bounds: (f32, f32),
    bins: usize,
) -> PyResult<Vec<PyLinearDensityBin>> {
    let positions = borrowed_coordinates(&positions)?;
    let weight_values = weights.as_slice()?;
    py.detach(move || {
        pdbiox::analysis::linear_density(
            positions,
            weight_values,
            pdbiox::analysis::LinearDensityOptions {
                axis: axis.into(),
                minimum: bounds.0,
                maximum: bounds.1,
                bins,
            },
        )
    })
    .map(|values| values.into_iter().map(PyLinearDensityBin::from).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn linear_density_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    weights: PyReadonlyArray1<'_, f64>,
    options: PyLinearDensityOptions,
) -> PyResult<Vec<PyLinearDensityBin>> {
    let positions = borrowed_coordinates(&positions)?;
    let weights = weights.as_slice()?;
    py.detach(move || {
        pdbiox::analysis::linear_density(
            positions,
            weights,
            pdbiox::analysis::LinearDensityOptions {
                axis: options.axis.into(),
                minimum: options.minimum,
                maximum: options.maximum,
                bins: options.bins,
            },
        )
    })
    .map(|values| values.into_iter().map(PyLinearDensityBin::from).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn density_map(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    weights: PyReadonlyArray1<'_, f64>,
    spec: &PyDensityGridSpec,
) -> PyResult<PyDensityGrid> {
    let positions = borrowed_coordinates(&positions)?;
    let weight_values = weights.as_slice()?;
    let spec = spec.0;
    py.detach(move || pdbiox::analysis::density_map(positions, weight_values, spec))
        .map(|value| PyDensityGrid {
            spec: value.spec,
            values: value.density,
            excluded_weight: value.excluded_weight,
        })
        .map_err(value_error)
}

impl From<PyCartesianAxis> for pdbiox::analysis::CartesianAxis {
    fn from(value: PyCartesianAxis) -> Self {
        match value {
            PyCartesianAxis::X => Self::X,
            PyCartesianAxis::Y => Self::Y,
            PyCartesianAxis::Z => Self::Z,
        }
    }
}

impl From<pdbiox::analysis::CartesianAxis> for PyCartesianAxis {
    fn from(value: pdbiox::analysis::CartesianAxis) -> Self {
        match value {
            pdbiox::analysis::CartesianAxis::X => Self::X,
            pdbiox::analysis::CartesianAxis::Y => Self::Y,
            pdbiox::analysis::CartesianAxis::Z => Self::Z,
        }
    }
}

impl From<pdbiox::analysis::LinearDensityBin> for PyLinearDensityBin {
    fn from(value: pdbiox::analysis::LinearDensityBin) -> Self {
        Self {
            lower: value.lower,
            upper: value.upper,
            weight: value.weight,
            density: value.density,
        }
    }
}

pub(crate) fn linear_density_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::LinearDensityBin>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        let values = values
            .into_iter()
            .map(|value| Py::new(py, PyLinearDensityBin::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, values)?.unbind().into_any())
    })
}

pub(crate) fn density_map_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<pdbiox::analysis::DensityGrid>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        Ok(Py::new(
            py,
            PyDensityGrid {
                spec: value.spec,
                values: value.density,
                excluded_weight: value.excluded_weight,
            },
        )?
        .into_any())
    })
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
