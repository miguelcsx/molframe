//! Python-owned result and mesh types for the reusable surface kernels.

use crate::geometry::coordinates;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "SurfacePoint", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfacePoint {
    #[pyo3(get)]
    atom: usize,
    #[pyo3(get)]
    position: [f32; 3],
    #[pyo3(get)]
    normal: [f32; 3],
}

#[pyclass(name = "ExcludedSurfacePoint", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyExcludedSurfacePoint {
    #[pyo3(get)]
    atom: usize,
    #[pyo3(get)]
    excluded: usize,
    #[pyo3(get)]
    position: [f32; 3],
    #[pyo3(get)]
    area: f64,
}

#[pyclass(name = "AtomContactArea", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomContactArea {
    #[pyo3(get)]
    first: usize,
    #[pyo3(get)]
    second: usize,
    #[pyo3(get)]
    area: f64,
}

#[pyclass(name = "AtomDepthOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomDepthOptions(pub(crate) pdbiox::surface::AtomDepthOptions);

#[pymethods]
impl PyAtomDepthOptions {
    #[new]
    fn new(cell_size: f64) -> Self {
        Self(pdbiox::surface::AtomDepthOptions { cell_size })
    }
}

#[pyclass(name = "BuriedSurface", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBuriedSurface {
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

#[pyclass(name = "MoleculeRole", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMoleculeRole {
    First,
    Second,
    Excluded,
}

#[pyclass(name = "SurfaceTriangle", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceTriangle {
    #[pyo3(get)]
    vertices: [[f32; 3]; 3],
    #[pyo3(get)]
    normal: [f32; 3],
    #[pyo3(get)]
    area: f64,
}

#[pyclass(name = "SolventExcludedSurface", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySolventExcludedSurface {
    #[pyo3(get)]
    triangles: Vec<PySurfaceTriangle>,
    #[pyo3(get)]
    area: f64,
    pub(crate) native: pdbiox::surface::SolventExcludedSurface,
}

#[pymethods]
impl PySolventExcludedSurface {
    fn indexed_mesh(&self) -> PyIndexedSurfaceMesh {
        PyIndexedSurfaceMesh::from_native(self.native.indexed_mesh())
    }
}

#[pyclass(name = "SurfaceFace", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceFace {
    #[pyo3(get)]
    indices: [u32; 3],
}

#[pymethods]
impl PySurfaceFace {
    #[new]
    fn new(indices: [u32; 3]) -> Self {
        Self { indices }
    }
}

#[pyclass(name = "MeshReport", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMeshReport {
    #[pyo3(get)]
    boundary_edges: u32,
    #[pyo3(get)]
    non_manifold_edges: u32,
    #[pyo3(get)]
    degenerate_faces: u32,
    #[pyo3(get)]
    components: u32,
}

#[pyclass(name = "IndexedSurfaceMesh", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyIndexedSurfaceMesh {
    pub(crate) native: pdbiox::surface::IndexedSurfaceMesh,
}

#[pymethods]
impl PyIndexedSurfaceMesh {
    #[new]
    fn new(vertices: PyReadonlyArray2<'_, f32>, faces: Vec<[u32; 3]>) -> PyResult<Self> {
        let vertices = coordinates(vertices)?;
        Ok(Self {
            native: pdbiox::surface::IndexedSurfaceMesh::new(
                vertices,
                faces
                    .into_iter()
                    .map(pdbiox::surface::SurfaceFace)
                    .collect(),
            ),
        })
    }

    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        array2(py, &self.native.vertices)
    }

    #[getter]
    fn faces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let faces = self
            .native
            .faces
            .iter()
            .map(|face| face.0)
            .collect::<Vec<_>>();
        array2(py, &faces)
    }

    #[getter]
    fn vertex_normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        array2(py, &self.native.vertex_normals)
    }

    #[getter]
    fn report(&self) -> PyMeshReport {
        self.native.report.into()
    }

    fn neighbours(&self, vertex: u32) -> Vec<u32> {
        self.native.neighbours(vertex)
    }
}

#[pyclass(name = "SurfaceDistances", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySurfaceDistances {
    #[pyo3(get)]
    pub(crate) source: u32,
    #[pyo3(get)]
    pub(crate) distances: Vec<f64>,
}

#[pyclass(name = "CurvatureQuality", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCurvatureQuality {
    Interior,
    Boundary,
    Degenerate,
}

#[pyclass(name = "SurfaceCurvature", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceCurvature {
    #[pyo3(get)]
    mean: f64,
    #[pyo3(get)]
    gaussian: f64,
    #[pyo3(get)]
    maximum: f64,
    #[pyo3(get)]
    minimum: f64,
    #[pyo3(get)]
    shape_index: f64,
    #[pyo3(get)]
    curvedness: f64,
    #[pyo3(get)]
    quality: PyCurvatureQuality,
}

fn array2<'py, T>(py: Python<'py>, values: &[[T; 3]]) -> PyResult<Bound<'py, PyArray2<T>>>
where
    T: numpy::Element + Copy,
{
    Array2::from_shape_vec(
        (values.len(), 3),
        values
            .iter()
            .flat_map(|value| value.iter().copied())
            .collect(),
    )
    .map(|value| value.into_pyarray(py))
    .map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

impl From<PyMoleculeRole> for pdbiox::surface::MoleculeRole {
    fn from(value: PyMoleculeRole) -> Self {
        match value {
            PyMoleculeRole::First => Self::First,
            PyMoleculeRole::Second => Self::Second,
            PyMoleculeRole::Excluded => Self::Excluded,
        }
    }
}

impl From<pdbiox::surface::SurfacePoint> for PySurfacePoint {
    fn from(value: pdbiox::surface::SurfacePoint) -> Self {
        Self {
            atom: value.atom,
            position: value.position,
            normal: value.normal,
        }
    }
}

impl From<pdbiox::surface::ExcludedSurfacePoint> for PyExcludedSurfacePoint {
    fn from(value: pdbiox::surface::ExcludedSurfacePoint) -> Self {
        Self {
            atom: value.atom,
            excluded: value.excluded,
            position: value.position,
            area: value.area,
        }
    }
}

impl From<pdbiox::surface::AtomContactArea> for PyAtomContactArea {
    fn from(value: pdbiox::surface::AtomContactArea) -> Self {
        Self {
            first: value.first,
            second: value.second,
            area: value.area,
        }
    }
}

impl From<pdbiox::surface::BuriedSurface> for PyBuriedSurface {
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

impl From<pdbiox::surface::SurfaceTriangle> for PySurfaceTriangle {
    fn from(value: pdbiox::surface::SurfaceTriangle) -> Self {
        Self {
            vertices: value.vertices,
            normal: value.normal,
            area: value.area,
        }
    }
}

impl From<pdbiox::surface::SolventExcludedSurface> for PySolventExcludedSurface {
    fn from(value: pdbiox::surface::SolventExcludedSurface) -> Self {
        Self {
            triangles: value.triangles.iter().copied().map(Into::into).collect(),
            area: value.area,
            native: value,
        }
    }
}

impl PyIndexedSurfaceMesh {
    pub(crate) fn from_native(native: pdbiox::surface::IndexedSurfaceMesh) -> Self {
        Self { native }
    }
}

impl From<pdbiox::surface::MeshReport> for PyMeshReport {
    fn from(value: pdbiox::surface::MeshReport) -> Self {
        Self {
            boundary_edges: value.boundary_edges,
            non_manifold_edges: value.non_manifold_edges,
            degenerate_faces: value.degenerate_faces,
            components: value.components,
        }
    }
}

impl From<pdbiox::surface::SurfaceDistances> for PySurfaceDistances {
    fn from(value: pdbiox::surface::SurfaceDistances) -> Self {
        Self {
            source: value.source,
            distances: value.distances.into_vec(),
        }
    }
}

impl From<pdbiox::surface::CurvatureQuality> for PyCurvatureQuality {
    fn from(value: pdbiox::surface::CurvatureQuality) -> Self {
        match value {
            pdbiox::surface::CurvatureQuality::Interior => Self::Interior,
            pdbiox::surface::CurvatureQuality::Boundary => Self::Boundary,
            pdbiox::surface::CurvatureQuality::Degenerate => Self::Degenerate,
        }
    }
}

impl From<pdbiox::surface::SurfaceCurvature> for PySurfaceCurvature {
    fn from(value: pdbiox::surface::SurfaceCurvature) -> Self {
        Self {
            mean: value.mean,
            gaussian: value.gaussian,
            maximum: value.maximum,
            minimum: value.minimum,
            shape_index: value.shape_index,
            curvedness: value.curvedness,
            quality: value.quality.into(),
        }
    }
}

pub(crate) fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySurfacePoint>()?;
    module.add_class::<PyExcludedSurfacePoint>()?;
    module.add_class::<PyAtomContactArea>()?;
    module.add_class::<PyAtomDepthOptions>()?;
    module.add_class::<PyBuriedSurface>()?;
    module.add_class::<PyCavity>()?;
    module.add_class::<PyMoleculeRole>()?;
    module.add_class::<PySurfaceTriangle>()?;
    module.add_class::<PySolventExcludedSurface>()?;
    module.add_class::<PySurfaceFace>()?;
    module.add_class::<PyMeshReport>()?;
    module.add_class::<PyIndexedSurfaceMesh>()?;
    module.add_class::<PySurfaceDistances>()?;
    module.add_class::<PyCurvatureQuality>()?;
    module.add_class::<PySurfaceCurvature>()?;
    Ok(())
}
