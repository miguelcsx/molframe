//! Borrowed numerical inputs to caller-installed analytical surface executables.
use molframe::surface::native_solver::{
    NativeSurfaceMesh, NativeSurfaceOptions, NativeSurfaceProvider, SurfaceAtoms,
    native_molecular_surface,
};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use std::{path::PathBuf, sync::Arc};

/// A caller-installed scientific executable and its actual surface controls.
#[derive(Clone, Debug)]
#[pyclass(
    name = "NativeSurfaceProvider",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PyNativeSurfaceProvider {
    executable: PathBuf,
    provider: NativeSurfaceProvider,
    options: NativeSurfaceOptions,
}

#[pymethods]
impl PyNativeSurfaceProvider {
    #[new]
    #[pyo3(signature = (executable, method, *, probe=1.4, density=3.0, all_components=true))]
    fn new(
        executable: PathBuf,
        method: &str,
        probe: f32,
        density: f32,
        all_components: bool,
    ) -> PyResult<Self> {
        let provider = match method {
            "surf" => NativeSurfaceProvider::Surf,
            "msms" => NativeSurfaceProvider::Msms,
            _ => return Err(crate::error::value("method must be surf or msms")),
        };
        if !probe.is_finite() || probe <= 0.0 || !density.is_finite() || density <= 0.0 {
            return Err(crate::error::value(
                "probe and density must be finite and positive",
            ));
        }
        Ok(Self {
            executable,
            provider,
            options: NativeSurfaceOptions {
                probe_radius: probe,
                density,
                all_components,
                ..NativeSurfaceOptions::default()
            },
        })
    }
}

/// Immutable provider output retaining its execution-memory reservation.
#[derive(Clone, Debug)]
#[pyclass(
    name = "MolecularSurface",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PyMolecularSurface {
    inner: Arc<NativeSurfaceMesh>,
}

fn rows3<T: numpy::Element>(
    py: Python<'_>,
    count: usize,
    values: Vec<T>,
) -> PyResult<Bound<'_, PyArray2<T>>> {
    let array = values.into_pyarray(py).reshape((count, 3))?;
    array.readwrite().make_nonwriteable();
    Ok(array)
}

#[pymethods]
impl PyMolecularSurface {
    /// Exactly zero-area provider triangles omitted without changing their boundary.
    #[getter]
    fn discarded_degenerate_faces(&self) -> u64 {
        self.inner.mesh.discarded_degenerate_faces
    }
    /// Generated surface coordinates, rather than source atom coordinates.
    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let column = &self.inner.mesh.geometry.vertices;
        rows3(py, column.len(), column.iter().flatten().copied().collect())
    }
    /// Provider's analytic unit normals.
    #[getter]
    fn normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let column = &self.inner.mesh.geometry.vertex_normals;
        rows3(py, column.len(), column.iter().flatten().copied().collect())
    }
    /// Triangle vertex indices.
    #[getter]
    fn faces<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let column = &self.inner.mesh.geometry.faces;
        rows3(
            py,
            column.len(),
            column.iter().flat_map(|face| face.0).collect(),
        )
    }
    /// Original input atom row at every generated vertex.
    #[getter]
    fn atom_rows<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let array = self.inner.mesh.atom_rows.to_vec().into_pyarray(py);
        array.readwrite().make_nonwriteable();
        array
    }
    /// Actual provider.
    #[getter]
    fn method(&self) -> &'static str {
        match self.inner.provider() {
            NativeSurfaceProvider::Surf => "surf",
            NativeSurfaceProvider::Msms => "msms",
        }
    }
    /// Actual rolling-probe radius, angstroms.
    #[getter]
    fn probe_radius(&self) -> f32 {
        self.inner.options().probe_radius
    }
    /// Actual MSMS vertices per square angstrom; unused by Surf.
    #[getter]
    fn density(&self) -> f32 {
        self.inner.options().density
    }
    /// Whether disconnected components were requested.
    #[getter]
    fn all_components(&self) -> bool {
        self.inner.options().all_components
    }
}

/// Computes an authentic Surf or MSMS result using the explicitly selected executable.
/// The execution context must authorize temporary storage. Ctrl-C cancels and reaps the child.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, rows, provider, *, context=None))]
fn native_surface(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    rows: &Bound<'_, PyArray1<u32>>,
    provider: &PyNativeSurfaceProvider,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyMolecularSurface> {
    let coordinates = coordinates.readonly();
    let radii = radii.readonly();
    let rows = rows.readonly();
    let positions = crate::bindings::coordinates(&coordinates)?;
    let radii = radii
        .as_slice()
        .map_err(|_| crate::error::value("radii must be C-contiguous"))?;
    let rows = rows
        .as_slice()
        .map_err(|_| crate::error::value("rows must be C-contiguous"))?;
    let result = crate::execution::run(py, context, |context| {
        native_molecular_surface(
            &provider.executable,
            provider.provider,
            SurfaceAtoms {
                positions,
                radii,
                rows,
            },
            provider.options,
            context,
        )
    })?
    .map_err(|error| provider_error(&error))?;
    Ok(PyMolecularSurface {
        inner: Arc::new(result),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNativeSurfaceProvider>()?;
    module.add_class::<PyMolecularSurface>()?;
    module.add_function(wrap_pyfunction!(native_surface, module)?)
}

fn provider_error(error: &molframe::surface::native_solver::NativeSurfaceError) -> PyErr {
    use molframe::surface::{SurfaceOutputError, native_solver::NativeSurfaceError};
    use molframe::{Code, Diagnostic};
    let code = match error {
        NativeSurfaceError::Cancelled(_)
        | NativeSurfaceError::Output(SurfaceOutputError::Cancelled(_)) => Code::E1904,
        NativeSurfaceError::Memory(_) => Code::E7001,
        NativeSurfaceError::Storage(_)
        | NativeSurfaceError::Timeout
        | NativeSurfaceError::Output(SurfaceOutputError::Limit(_)) => Code::E1903,
        NativeSurfaceError::Io(_)
        | NativeSurfaceError::Failed { .. }
        | NativeSurfaceError::Output(SurfaceOutputError::Io(_)) => Code::E7901,
        NativeSurfaceError::Invalid(_)
        | NativeSurfaceError::Output(SurfaceOutputError::Invalid(_)) => {
            return crate::error::value(error.to_string());
        }
    };
    crate::error::kernel(Diagnostic::new(code).with_message(error.to_string()))
}
