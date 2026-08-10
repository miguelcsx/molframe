//! Governed intrinsic workflows backed by Rust `Analysis<T>` results.

use super::{
    PyCurvature, PySurfaceGridOptions, PySurfaceMesh, array2, coordinate_frames_owned,
    diffusion_result, pca_result, scalar_values, triples_owned, value_error,
};
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::PyAnalysisPolicy;
use numpy::PyUntypedArrayMethods;
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pdbiox::PeriodicAngle;
use pdbiox::surface::{SurfaceWorkflowOptions, SurfaceWorkflowResult, governed_surface_geometry};
use pdbiox::traj::{
    CartesianFit, EnsembleDistanceMatrix, analyse_cartesian_pca, analyse_diffusion_map,
    analyse_dihedral_pca,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[derive(Clone, Debug)]
enum OwnedCartesianFit {
    None,
    Reference(Box<[[f32; 3]]>),
    IterativeMean {
        max_iterations: usize,
        tolerance: f64,
    },
}

#[pyclass(name = "CartesianFit", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCartesianFit(OwnedCartesianFit);

#[pymethods]
impl PyCartesianFit {
    #[staticmethod]
    fn none() -> Self {
        Self(OwnedCartesianFit::None)
    }

    #[staticmethod]
    fn reference(coordinates: PyReadonlyArray2<'_, f32>) -> PyResult<Self> {
        triples_owned(coordinates).map(|value| Self(OwnedCartesianFit::Reference(value)))
    }

    #[staticmethod]
    fn iterative_mean(max_iterations: usize, tolerance: f64) -> Self {
        Self(OwnedCartesianFit::IterativeMean {
            max_iterations,
            tolerance,
        })
    }
}

impl PyCartesianFit {
    pub(super) fn as_rust(&self) -> CartesianFit<'_> {
        match &self.0 {
            OwnedCartesianFit::None => CartesianFit::None,
            OwnedCartesianFit::Reference(reference) => CartesianFit::Reference(reference),
            OwnedCartesianFit::IterativeMean {
                max_iterations,
                tolerance,
            } => CartesianFit::IterativeMean {
                max_iterations: *max_iterations,
                tolerance: *tolerance,
            },
        }
    }
}

#[pyclass(name = "SurfaceWorkflowResult", frozen)]
pub(crate) struct PySurfaceWorkflowResult {
    #[pyo3(get)]
    mesh: Py<PySurfaceMesh>,
    #[pyo3(get)]
    curvature: Vec<PyCurvature>,
    #[pyo3(get)]
    geodesic_source: Option<u32>,
    geodesic_distances: Option<Py<PyArray1<f64>>>,
}

#[pyclass(name = "SurfaceWorkflowOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySurfaceWorkflowOptions {
    selection: String,
    radii_set: String,
    probe: f32,
    grid: PySurfaceGridOptions,
    source_vertex: Option<u32>,
}

#[pymethods]
impl PySurfaceWorkflowOptions {
    #[new]
    #[pyo3(signature = (selection, radii_set, probe, grid, source_vertex=None))]
    fn new(
        selection: String,
        radii_set: String,
        probe: f32,
        grid: PySurfaceGridOptions,
        source_vertex: Option<u32>,
    ) -> Self {
        Self {
            selection,
            radii_set,
            probe,
            grid,
            source_vertex,
        }
    }
}

#[pymethods]
impl PySurfaceWorkflowResult {
    #[getter]
    fn geodesic_distances(&self, py: Python<'_>) -> Option<Py<PyArray1<f64>>> {
        self.geodesic_distances
            .as_ref()
            .map(|value| value.clone_ref(py))
    }
}

#[pyfunction]
#[pyo3(signature = (positions, radii, options, policy))]
pub(crate) fn analyse_surface_geometry(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    options: &PySurfaceWorkflowOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let positions = triples_owned(positions)?;
    let radii = scalar_values(radii);
    let options = options.clone();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            governed_surface_geometry(
                &positions,
                &radii,
                SurfaceWorkflowOptions {
                    selection: &options.selection,
                    radii_set: &options.radii_set,
                    probe: options.probe,
                    grid: options.grid.0,
                    source_vertex: options.source_vertex,
                },
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, surface_workflow_value)
}

#[pyfunction]
#[pyo3(signature = (frames, *, fit, components, memory_limit, policy))]
pub(crate) fn analyse_pca(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    fit: &PyCartesianFit,
    components: usize,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames = coordinate_frames_owned(frames);
    let fit = fit.clone();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| analyse_cartesian_pca(&frames, fit.as_rust(), components, memory_limit, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, result, |py, value| {
        Py::new(py, pca_result(py, value).map_err(value_error)?).map(Py::into_any)
    })
}

#[pyfunction]
#[pyo3(signature = (angles, *, torsion_set, components, memory_limit, policy))]
pub(crate) fn analyse_torsion_pca(
    py: Python<'_>,
    angles: PyReadonlyArray2<'_, f64>,
    torsion_set: &str,
    components: usize,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let periodic = periodic_angles(angles)?;
    let torsion_set = torsion_set.to_owned();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| analyse_dihedral_pca(&periodic, &torsion_set, components, memory_limit, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, result, |py, value| {
        Py::new(py, pca_result(py, value).map_err(value_error)?).map(Py::into_any)
    })
}

#[pyfunction]
#[pyo3(signature = (distances, *, metric, epsilon, time, dimensions, policy))]
pub(crate) fn analyse_diffusion(
    py: Python<'_>,
    distances: PyReadonlyArray2<'_, f64>,
    metric: &str,
    epsilon: f64,
    time: u32,
    dimensions: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let matrix = distance_matrix_owned(distances)?;
    let metric = metric.to_owned();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| analyse_diffusion_map(&matrix, &metric, epsilon, time, dimensions, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, result, |py, value| {
        Py::new(py, diffusion_result(py, value).map_err(value_error)?).map(Py::into_any)
    })
}

fn periodic_angles(values: PyReadonlyArray2<'_, f64>) -> PyResult<Vec<Vec<PeriodicAngle>>> {
    let result = values
        .as_array()
        .outer_iter()
        .map(|row| {
            row.iter()
                .copied()
                .map(PeriodicAngle::from_radians)
                .collect::<Result<Vec<_>, _>>()
                .map_err(value_error)
        })
        .collect();
    drop(values);
    result
}

fn distance_matrix_owned(values: PyReadonlyArray2<'_, f64>) -> PyResult<EnsembleDistanceMatrix> {
    let shape = values.shape();
    if shape[0] != shape[1] {
        return Err(PyValueError::new_err("distance matrix must be square"));
    }
    let matrix = EnsembleDistanceMatrix {
        size: shape[0],
        values: values.as_array().iter().copied().collect(),
    };
    drop(values);
    Ok(matrix)
}

fn surface_workflow_value(py: Python<'_>, value: SurfaceWorkflowResult) -> PyResult<Py<PyAny>> {
    let curvature = value
        .curvature
        .into_iter()
        .map(|item| {
            (
                item.mean,
                item.gaussian,
                item.maximum,
                item.minimum,
                item.shape_index,
                item.curvedness,
                format!("{:?}", item.quality),
            )
        })
        .collect();
    let vertices = array2(py, &value.mesh.vertices).map_err(value_error)?;
    let faces = array2(
        py,
        &value
            .mesh
            .faces
            .iter()
            .map(|face| face.0)
            .collect::<Vec<_>>(),
    )
    .map_err(value_error)?;
    let mesh = Py::new(
        py,
        PySurfaceMesh {
            mesh: value.mesh,
            vertices,
            faces,
        },
    )?;
    let (geodesic_source, geodesic_distances) = value.geodesics.map_or((None, None), |distances| {
        (
            Some(distances.source),
            Some(
                Array1::from_vec(distances.distances.into_vec())
                    .into_pyarray(py)
                    .unbind(),
            ),
        )
    });
    Py::new(
        py,
        PySurfaceWorkflowResult {
            mesh,
            curvature,
            geodesic_source,
            geodesic_distances,
        },
    )
    .map(Py::into_any)
}
