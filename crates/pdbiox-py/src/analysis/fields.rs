//! Vectorized density and surface kernels with caller-owned scientific policy.

use crate::geometry::coordinates;
use crate::intrinsic::PySurfaceGridOptions;
use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray1, PyArray3, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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
    fn density<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        Array3::from_shape_vec(self.spec.shape, self.values.clone())
            .map(|array| array.into_pyarray(py))
            .map_err(value_error)
    }
}

#[pyclass(name = "SurfaceAreas", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySurfaceAreas {
    #[pyo3(get)]
    first_alone: f64,
    #[pyo3(get)]
    second_alone: f64,
    #[pyo3(get)]
    together: f64,
    #[pyo3(get)]
    buried: f64,
}

#[pyclass(name = "Cavity", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCavity {
    #[pyo3(get)]
    volume: f64,
    #[pyo3(get)]
    representative: [f32; 3],
    #[pyo3(get)]
    cells: usize,
}

#[pyclass(name = "AtomDepthOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomDepthOptions(pdbiox::surface::AtomDepthOptions);

#[pymethods]
impl PyAtomDepthOptions {
    #[new]
    fn new(cell_size: f64) -> Self {
        Self(pdbiox::surface::AtomDepthOptions { cell_size })
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
    let positions = coordinates(positions)?;
    let weight_values = weights.as_slice()?.to_vec();
    drop(weights);
    py.detach(move || {
        pdbiox::analysis::linear_density(
            &positions,
            &weight_values,
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
pub(crate) fn density_map(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    weights: PyReadonlyArray1<'_, f64>,
    spec: &PyDensityGridSpec,
) -> PyResult<PyDensityGrid> {
    let positions = coordinates(positions)?;
    let weight_values = weights.as_slice()?.to_vec();
    drop(weights);
    let spec = spec.0;
    py.detach(move || pdbiox::analysis::density_map(&positions, &weight_values, spec))
        .map(|value| PyDensityGrid {
            spec: value.spec,
            values: value.density,
            excluded_weight: value.excluded_weight,
        })
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn solvent_accessible_surface(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe_radius: f32,
    sample_points: u16,
) -> PyResult<Vec<f64>> {
    let positions = coordinates(positions)?;
    let radius_values = radii.as_slice()?.to_vec();
    drop(radii);
    py.detach(move || {
        pdbiox::surface::shrake_rupley(&positions, &radius_values, probe_radius, sample_points)
    })
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn buried_surface(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    first: PyReadonlyArray1<'_, bool>,
    probe_radius: f32,
    sample_points: u16,
) -> PyResult<PySurfaceAreas> {
    let positions = coordinates(positions)?;
    let radius_values = radii.as_slice()?.to_vec();
    drop(radii);
    let first_values = first.as_slice()?.to_vec();
    drop(first);
    py.detach(move || {
        pdbiox::surface::buried_surface(
            &positions,
            &radius_values,
            probe_radius,
            sample_points,
            &first_values,
        )
    })
    .map(PySurfaceAreas::from)
    .map_err(value_error)
}

/// Compute all atom depths through one native Rust kernel invocation.
#[pyfunction]
pub(crate) fn atom_depths<'py>(
    py: Python<'py>,
    atoms: PyReadonlyArray2<'_, f32>,
    surface: PyReadonlyArray2<'_, f32>,
    options: PyAtomDepthOptions,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let atoms = coordinates(atoms)?;
    let surface = coordinates(surface)?;
    let options = options.0;
    py.detach(move || pdbiox::surface::atom_depths(&atoms, &surface, options))
        .map(|values| values.into_pyarray(py))
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn cavities(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii_array: PyReadonlyArray1<'_, f32>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<Vec<PyCavity>> {
    let positions = coordinates(positions)?;
    let radii = radii_array.as_slice()?.to_vec();
    drop(radii_array);
    py.detach(move || pdbiox::surface::cavities_with_options(&positions, &radii, probe, options.0))
        .map(|values| values.into_iter().map(Into::into).collect())
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

impl From<pdbiox::surface::BuriedSurface> for PySurfaceAreas {
    fn from(value: pdbiox::surface::BuriedSurface) -> Self {
        Self {
            first_alone: value.first_alone,
            second_alone: value.second_alone,
            together: value.together,
            buried: value.buried,
        }
    }
}

impl From<pdbiox::surface::Cavity> for PyCavity {
    fn from(value: pdbiox::surface::Cavity) -> Self {
        Self {
            volume: value.volume,
            representative: value.representative,
            cells: value.cells,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
