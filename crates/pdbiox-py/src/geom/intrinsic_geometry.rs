//! Thin Python handles over intrinsic-geometry Rust results.

#[path = "intrinsic_geometry/arrays.rs"]
mod arrays;
#[path = "intrinsic_geometry/ensemble.rs"]
mod ensemble;
#[path = "intrinsic_geometry/governed.rs"]
mod governed;

pub(crate) use arrays::{
    array2, borrowed_frame_input, borrowed_frames, diffusion_result, pca_result, scalar_values,
    triples_owned,
};

use crate::geometry::PyRotationOptions;
pub(crate) use ensemble::*;
pub(crate) use governed::{
    PyCartesianFit, PySurfaceWorkflowOptions, PySurfaceWorkflowResult, analyse_diffusion,
    analyse_pca, analyse_surface_geometry, analyse_torsion_pca,
};
use numpy::{
    PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods,
};
use pdbiox::surface::{
    IndexedSurfaceMesh, SurfaceGridOptions, solvent_excluded_surface_with_options,
    surface_curvatures,
};
use pdbiox::traj::{
    EnsembleDistanceMatrix, cartesian_pca_view as rust_pca_view, diffusion_map as rust_diffusion,
    dihedral_pca,
};
use pdbiox::{PeriodicAngle, Rotation3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type PyCurvature = (f64, f64, f64, f64, f64, f64, String);

#[pyclass(name = "SurfaceGridOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceGridOptions(pub(crate) SurfaceGridOptions);

#[pymethods]
impl PySurfaceGridOptions {
    #[new]
    fn new(resolution: f32, max_cells: usize) -> Self {
        Self(SurfaceGridOptions {
            resolution,
            max_cells,
        })
    }

    #[staticmethod]
    pub(crate) fn standard(resolution: f32) -> Self {
        Self(SurfaceGridOptions::standard(resolution))
    }

    #[classattr]
    const STANDARD_MAX_CELLS: usize = SurfaceGridOptions::STANDARD_MAX_CELLS;
}

#[pyclass(name = "PeriodicAngle", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPeriodicAngle {
    value: PeriodicAngle,
}

#[pymethods]
impl PyPeriodicAngle {
    #[new]
    #[pyo3(signature = (radians))]
    fn new(radians: f64) -> PyResult<Self> {
        PeriodicAngle::from_radians(radians)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn from_degrees(degrees: f64) -> PyResult<Self> {
        PeriodicAngle::from_degrees(degrees)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn from_radians(radians: f64) -> PyResult<Self> {
        Self::new(radians)
    }

    #[getter]
    fn radians(&self) -> f64 {
        self.value.radians()
    }

    #[getter]
    fn degrees(&self) -> f64 {
        self.value.degrees()
    }

    fn delta_to(&self, other: &Self) -> f64 {
        self.value.signed_delta(other.value)
    }

    fn signed_delta(&self, other: &Self) -> f64 {
        self.value.signed_delta(other.value)
    }

    fn distance_to(&self, other: &Self) -> f64 {
        self.value.distance(other.value)
    }

    fn distance(&self, other: &Self) -> f64 {
        self.value.distance(other.value)
    }
}

#[pyclass(name = "Rotation3", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRotation3 {
    pub(crate) value: Rotation3,
}

impl PyRotation3 {
    pub(crate) const fn from_rust(value: Rotation3) -> Self {
        Self { value }
    }
}

#[pymethods]
impl PyRotation3 {
    #[classattr]
    const IDENTITY: Self = Self {
        value: Rotation3::IDENTITY,
    };

    #[new]
    fn new(matrix: [[f64; 3]; 3]) -> PyResult<Self> {
        Rotation3::from_matrix(matrix, 1.0e-10)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn exponential(rotation_vector: [f64; 3]) -> PyResult<Self> {
        Rotation3::exp(rotation_vector)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn identity() -> Self {
        Self {
            value: Rotation3::IDENTITY,
        }
    }

    #[staticmethod]
    fn from_matrix(matrix: [[f64; 3]; 3], tolerance: f64) -> PyResult<Self> {
        Rotation3::from_matrix(matrix, tolerance)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn exponential_with_options(
        rotation_vector: [f64; 3],
        options: &PyRotationOptions,
    ) -> PyResult<Self> {
        Rotation3::exp_with_options(rotation_vector, options.0)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    #[staticmethod]
    fn exp(rotation_vector: [f64; 3]) -> PyResult<Self> {
        Self::exponential(rotation_vector)
    }

    #[staticmethod]
    fn exp_with_options(rotation_vector: [f64; 3], options: &PyRotationOptions) -> PyResult<Self> {
        Self::exponential_with_options(rotation_vector, options)
    }

    #[getter]
    fn matrix(&self) -> [[f64; 3]; 3] {
        self.value.matrix()
    }

    fn logarithm(&self) -> [f64; 3] {
        self.value.log()
    }

    fn log(&self) -> [f64; 3] {
        self.value.log()
    }

    fn log_with_options(&self, options: &PyRotationOptions) -> PyResult<[f64; 3]> {
        self.value.log_with_options(options.0).map_err(value_error)
    }

    fn apply(&self, vector: [f64; 3]) -> [f64; 3] {
        self.value.apply(vector)
    }

    fn inverse(&self) -> Self {
        Self {
            value: self.value.inverse(),
        }
    }

    fn then(&self, next: &Self) -> Self {
        Self {
            value: self.value.then(next.value),
        }
    }

    fn distance_to(&self, other: &Self) -> f64 {
        self.value.distance(other.value)
    }

    fn interpolate(&self, other: &Self, fraction: f64) -> PyResult<Self> {
        self.value
            .interpolate(other.value, fraction)
            .map(|value| Self { value })
            .map_err(value_error)
    }

    fn interpolate_with_options(
        &self,
        other: &Self,
        fraction: f64,
        options: &PyRotationOptions,
    ) -> PyResult<Self> {
        self.value
            .interpolate_with_options(other.value, fraction, options.0)
            .map(|value| Self { value })
            .map_err(value_error)
    }
}

#[pyclass(name = "SurfaceMesh", frozen)]
pub(crate) struct PySurfaceMesh {
    mesh: IndexedSurfaceMesh,
    vertices: Py<PyArray2<f32>>,
    faces: Py<PyArray2<u32>>,
}

#[pymethods]
impl PySurfaceMesh {
    #[getter]
    fn vertices(&self, py: Python<'_>) -> Py<PyArray2<f32>> {
        self.vertices.clone_ref(py)
    }

    #[getter]
    fn faces(&self, py: Python<'_>) -> Py<PyArray2<u32>> {
        self.faces.clone_ref(py)
    }

    fn curvature(&self) -> PyResult<Vec<PyCurvature>> {
        surface_curvatures(&self.mesh)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        (
                            value.mean,
                            value.gaussian,
                            value.maximum,
                            value.minimum,
                            value.shape_index,
                            value.curvedness,
                            format!("{:?}", value.quality),
                        )
                    })
                    .collect()
            })
            .map_err(value_error)
    }
}

#[pyfunction]
pub(crate) fn surface_mesh(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<PySurfaceMesh> {
    let positions = triples_owned(positions)?;
    let radii = scalar_values(radii);
    let surface = py
        .detach(move || solvent_excluded_surface_with_options(&positions, &radii, probe, options.0))
        .map_err(value_error)?;
    let mesh = surface.indexed_mesh();
    let vertices = array2(py, &mesh.vertices).map_err(value_error)?;
    let faces = array2(
        py,
        &mesh.faces.iter().map(|face| face.0).collect::<Vec<_>>(),
    )
    .map_err(value_error)?;
    Ok(PySurfaceMesh {
        mesh,
        vertices,
        faces,
    })
}

#[pyclass(name = "PcaResult", frozen)]
pub(crate) struct PyPcaResult {
    mean: Py<PyArray1<f64>>,
    eigenvalues: Py<PyArray1<f64>>,
    components: Py<PyArray2<f64>>,
    projections: Py<PyArray2<f64>>,
}

#[pymethods]
impl PyPcaResult {
    #[getter]
    fn mean(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.mean.clone_ref(py)
    }

    #[getter]
    fn eigenvalues(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.eigenvalues.clone_ref(py)
    }

    #[getter]
    fn components(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.components.clone_ref(py)
    }

    #[getter]
    fn projections(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.projections.clone_ref(py)
    }
}

#[pyfunction]
#[pyo3(signature = (frames, fit, components=3, memory_limit=536_870_912))]
pub(crate) fn cartesian_pca(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    fit: &PyCartesianFit,
    components: usize,
    memory_limit: usize,
) -> PyResult<PyPcaResult> {
    let frames = borrowed_frames(&frames)?;
    let fit = fit.clone();
    let result = py
        .detach(|| rust_pca_view(frames, fit.as_rust(), components, memory_limit))
        .map_err(value_error)?;
    pca_result(py, result).map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (angles, components=3, memory_limit=536_870_912))]
pub(crate) fn torsion_pca(
    py: Python<'_>,
    angles: PyReadonlyArray2<'_, f64>,
    components: usize,
    memory_limit: usize,
) -> PyResult<PyPcaResult> {
    let periodic: Result<Vec<Vec<_>>, _> = angles
        .as_array()
        .outer_iter()
        .map(|row| {
            row.iter()
                .copied()
                .map(PeriodicAngle::from_radians)
                .collect()
        })
        .collect();
    drop(angles);
    let periodic = periodic.map_err(value_error)?;
    let result = py
        .detach(|| dihedral_pca(&periodic, components, memory_limit))
        .map_err(value_error)?;
    pca_result(py, result).map_err(value_error)
}

#[pyclass(name = "DiffusionMap", frozen)]
pub(crate) struct PyDiffusionMap {
    eigenvalues: Py<PyArray1<f64>>,
    coordinates: Py<PyArray2<f64>>,
    graph_components: Py<PyArray1<usize>>,
}

#[pymethods]
impl PyDiffusionMap {
    #[getter]
    fn eigenvalues(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.eigenvalues.clone_ref(py)
    }

    #[getter]
    fn coordinates(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.coordinates.clone_ref(py)
    }

    #[getter]
    fn graph_components(&self, py: Python<'_>) -> Py<PyArray1<usize>> {
        self.graph_components.clone_ref(py)
    }
}

#[pyfunction]
#[pyo3(signature = (distances, epsilon, time=1, dimensions=3))]
pub(crate) fn diffusion_map(
    py: Python<'_>,
    distances: PyReadonlyArray2<'_, f64>,
    epsilon: f64,
    time: u32,
    dimensions: usize,
) -> PyResult<PyDiffusionMap> {
    let shape = distances.shape();
    if shape[0] != shape[1] {
        return Err(PyValueError::new_err("distance matrix must be square"));
    }
    let matrix = EnsembleDistanceMatrix {
        size: shape[0],
        values: distances.as_array().iter().copied().collect(),
    };
    drop(distances);
    let result = py
        .detach(|| rust_diffusion(&matrix, epsilon, time, dimensions))
        .map_err(value_error)?;
    diffusion_result(py, result).map_err(value_error)
}

fn value_error(error: impl std::fmt::Debug) -> PyErr {
    PyValueError::new_err(format!("{error:?}"))
}
