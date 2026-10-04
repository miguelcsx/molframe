//! Triangle meshes of molecular surfaces: construction, components, curvature and
//! distances along the surface.

use molframe::surface::{
    CurvatureQuality, IndexedSurfaceMesh, SurfaceComponentFilter, SurfaceFace,
    edge_geodesic_distances, filter_surface_components, surface_components, surface_curvatures,
    surface_patch, write_obj,
};
use molframe::{Code, Diagnostic};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

/// `(rows, 3)` array from flat triples.
fn rows3<T: numpy::Element>(
    py: Python<'_>,
    rows: usize,
    flat: Vec<T>,
) -> PyResult<Bound<'_, PyArray2<T>>> {
    flat.into_pyarray(py).reshape((rows, 3))
}

fn frozen<T: numpy::Element>(array: Bound<'_, PyArray2<T>>) -> Bound<'_, PyArray2<T>> {
    array.readwrite().make_nonwriteable();
    array
}

/// An indexed triangle mesh: shared vertices and faces of three vertex indices.
///
/// The mesh reports, rather than repairs, what is wrong with it: boundary edges,
/// edges shared by more than two faces, degenerate faces (including a face that names
/// a vertex the mesh does not have). The distance and curvature operations refuse a mesh
/// that is not manifold.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Mesh",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PyMesh {
    pub(crate) inner: Arc<IndexedSurfaceMesh>,
}

impl PyMesh {
    pub(crate) fn new(mesh: IndexedSurfaceMesh) -> Self {
        Self {
            inner: Arc::new(mesh),
        }
    }
}

#[pymethods]
impl PyMesh {
    #[new]
    fn construct(
        vertices: &Bound<'_, PyArray2<f32>>,
        faces: &Bound<'_, PyArray2<u32>>,
    ) -> PyResult<Self> {
        let (vertices, faces) = (vertices.readonly(), faces.readonly());
        if vertices.as_array().ncols() != 3 || faces.as_array().ncols() != 3 {
            return Err(crate::error::value(
                "vertices and faces must each have three columns",
            ));
        }
        let vertices = vertices
            .as_array()
            .rows()
            .into_iter()
            .map(|row| [row[0], row[1], row[2]])
            .collect();
        let faces = faces
            .as_array()
            .rows()
            .into_iter()
            .map(|row| SurfaceFace([row[0], row[1], row[2]]))
            .collect();
        Ok(Self::new(IndexedSurfaceMesh::new(vertices, faces)))
    }

    /// Vertex positions, `(n, 3)`.
    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat = self.inner.vertices.iter().flatten().copied().collect();
        rows3(py, self.inner.vertices.len(), flat).map(frozen)
    }

    /// Triangles as vertex indices, `(m, 3)`.
    #[getter]
    fn faces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let flat = self.inner.faces.iter().flat_map(|face| face.0).collect();
        rows3(py, self.inner.faces.len(), flat).map(frozen)
    }

    /// Unit vertex normals, `(n, 3)`.
    #[getter]
    fn vertex_normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat = self
            .inner
            .vertex_normals
            .iter()
            .flatten()
            .copied()
            .collect();
        rows3(py, self.inner.vertex_normals.len(), flat).map(frozen)
    }

    /// Total area of the faces, in square ångström when the vertices are in ångström.
    #[getter]
    fn area(&self) -> f64 {
        self.inner.area()
    }

    /// Whether every edge has one or two faces and no face is degenerate.
    #[getter]
    fn is_manifold(&self) -> bool {
        self.inner.report.is_manifold()
    }

    /// Edges with exactly one face.
    #[getter]
    fn boundary_edges(&self) -> u32 {
        self.inner.report.boundary_edges
    }

    /// Edges shared by more than two faces.
    #[getter]
    fn non_manifold_edges(&self) -> u32 {
        self.inner.report.non_manifold_edges
    }

    /// Faces with a repeated or missing vertex, or no area.
    #[getter]
    fn degenerate_faces(&self) -> u32 {
        self.inner.report.degenerate_faces
    }

    /// Vertex-connected components.
    #[getter]
    fn component_count(&self) -> u32 {
        self.inner.report.components
    }

    /// Face-connected components as `(faces, area)`, in stable source-face order.
    fn components<'py>(&self, py: Python<'py>) -> Vec<(Bound<'py, PyArray1<u32>>, f64)> {
        surface_components(&self.inner)
            .into_iter()
            .map(|component| (component.faces.into_pyarray(py), component.area))
            .collect()
    }

    /// The mesh keeping the largest components of at least `minimum_area`, at most
    /// `maximum_components` of them.
    #[pyo3(signature = (*, minimum_area=0.0, maximum_components=None))]
    fn filter_components(
        &self,
        minimum_area: f64,
        maximum_components: Option<usize>,
    ) -> PyResult<Self> {
        filter_surface_components(
            &self.inner,
            SurfaceComponentFilter {
                minimum_area,
                maximum_components,
            },
        )
        .map(Self::new)
        .map_err(crate::error::kernel)
    }

    /// Curvature descriptors at every vertex; a non-manifold mesh is refused.
    fn curvatures(&self) -> PyResult<PyCurvatures> {
        let found = surface_curvatures(&self.inner).map_err(crate::error::kernel)?;
        Ok(PyCurvatures { found })
    }

    /// Distance from `source` to every vertex along the mesh edges.
    ///
    /// A graph approximation of the geodesic distance, not the exact one; infinity marks
    /// another component.
    fn geodesic_distances<'py>(
        &self,
        py: Python<'py>,
        source: u32,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let found = edge_geodesic_distances(&self.inner, source).map_err(crate::error::kernel)?;
        Ok(found.distances.into_vec().into_pyarray(py))
    }

    /// Vertices within `radius` of `source` along the mesh edges.
    fn patch<'py>(
        &self,
        py: Python<'py>,
        source: u32,
        radius: f64,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let found = surface_patch(&self.inner, source, radius).map_err(crate::error::kernel)?;
        Ok(found.into_vec().into_pyarray(py))
    }

    /// Writes Wavefront OBJ with vertex normals; an existing file is never replaced.
    #[allow(clippy::needless_pass_by_value)]
    fn write_obj(&self, path: PathBuf) -> PyResult<()> {
        write_obj(&path, &self.inner).map_err(|error| {
            crate::error::kernel(Diagnostic::new(Code::E7901).with_message(error.to_string()))
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Mesh(vertices={}, faces={}, manifold={})",
            self.inner.vertices.len(),
            self.inner.faces.len(),
            self.inner.report.is_manifold()
        )
    }
}

/// Curvature at every vertex of a mesh.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Curvatures",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PyCurvatures {
    found: Box<[molframe::surface::SurfaceCurvature]>,
}

impl PyCurvatures {
    fn column<'py>(
        &self,
        py: Python<'py>,
        pick: impl Fn(&molframe::surface::SurfaceCurvature) -> f64,
    ) -> Bound<'py, PyArray1<f64>> {
        self.found
            .iter()
            .map(pick)
            .collect::<Vec<_>>()
            .into_pyarray(py)
    }
}

#[pymethods]
impl PyCurvatures {
    /// Mean curvature.
    #[getter]
    fn mean<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.mean)
    }

    /// Gaussian curvature.
    #[getter]
    fn gaussian<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.gaussian)
    }

    /// Largest principal curvature.
    #[getter]
    fn maximum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.maximum)
    }

    /// Smallest principal curvature.
    #[getter]
    fn minimum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.minimum)
    }

    /// Koenderink shape index.
    #[getter]
    fn shape_index<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.shape_index)
    }

    /// Curvedness.
    #[getter]
    fn curvedness<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, |curvature| curvature.curvedness)
    }

    /// How far each estimate can be trusted: `interior`, `boundary` or `degenerate`.
    #[getter]
    fn quality(&self) -> Vec<&'static str> {
        self.found
            .iter()
            .map(|curvature| match curvature.quality {
                CurvatureQuality::Interior => "interior",
                CurvatureQuality::Boundary => "boundary",
                CurvatureQuality::Degenerate => "degenerate",
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.found.len()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMesh>()?;
    module.add_class::<PyCurvatures>()
}
