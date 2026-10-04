//! Mechanical adapters for solvent-accessible area and cavity detection.

use molframe::surface::{AtomDepthOptions, MoleculeRole, SurfaceGridOptions};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, ToPyArray};
use pyo3::prelude::*;

fn radii<'a>(array: &'a PyReadonlyArray1<'_, f32>) -> PyResult<&'a [f32]> {
    array.as_slice().map_err(|_| {
        crate::error::value("radii must be C-contiguous; call numpy.ascontiguousarray")
    })
}

/// Per-atom solvent-accessible area in square ångström (Shrake–Rupley).
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, points=960, context=None))]
fn sasa<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    radii: &Bound<'py, PyArray1<f32>>,
    probe: f32,
    points: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let areas = crate::execution::run(py, context, |context| {
        molframe::surface::shrake_rupley(positions, radii, probe, points, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(areas.to_pyarray(py))
}

/// Per-atom solvent-accessible area in square ångström (Lee–Richards slices).
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, slices=20, context=None))]
fn lee_richards<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    radii: &Bound<'py, PyArray1<f32>>,
    probe: f32,
    slices: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let areas = crate::execution::run(py, context, |context| {
        molframe::surface::lee_richards(positions, radii, probe, slices, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(areas.to_pyarray(py))
}

/// Enclosed cavities, largest first, as `(volume, (x, y, z), cells)` rows.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, resolution=0.5))]
fn cavities(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    probe: f32,
    resolution: f32,
) -> PyResult<Vec<(f64, [f32; 3], usize)>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let found = py
        .detach(|| molframe::surface::cavities(positions, radii, probe, resolution))
        .map_err(crate::error::kernel)?;
    Ok(found
        .into_iter()
        .map(|cavity| (cavity.volume, cavity.representative, cavity.cells))
        .collect())
}

/// How a surface calculation is bounded on its grid: the cell edge and the ceilings.
fn grid(
    resolution: f32,
    max_cells: Option<usize>,
    max_workspace_bytes: Option<usize>,
) -> SurfaceGridOptions {
    let standard = SurfaceGridOptions::standard(resolution);
    SurfaceGridOptions {
        resolution,
        max_cells: match max_cells {
            Some(cells) => cells,
            None => standard.max_cells,
        },
        max_workspace_bytes: match max_workspace_bytes {
            Some(bytes) => bytes,
            None => standard.max_workspace_bytes,
        },
    }
}

/// The surface points of one calculation: where, which way is out, on which atom.
#[derive(Clone, Debug)]
#[pyclass(
    name = "SurfacePoints",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PySurfacePoints {
    points: Vec<molframe::surface::SurfacePoint>,
}

#[pymethods]
impl PySurfacePoints {
    /// Point positions, `(n, 3)`.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat: Vec<f32> = self
            .points
            .iter()
            .flat_map(|point| point.position)
            .collect();
        flat.into_pyarray(py).reshape((self.points.len(), 3))
    }

    /// Outward unit normals, `(n, 3)`.
    #[getter]
    fn normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat: Vec<f32> = self.points.iter().flat_map(|point| point.normal).collect();
        flat.into_pyarray(py).reshape((self.points.len(), 3))
    }

    /// The atom each point sits on.
    #[getter]
    fn atoms<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        let atoms: Vec<u64> = self.points.iter().map(|point| point.atom as u64).collect();
        atoms.into_pyarray(py)
    }

    fn __len__(&self) -> usize {
        self.points.len()
    }
}

/// The accessible surface as points with outward normals, `samples` test directions per atom.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, samples=960, context=None))]
fn surface_points(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    probe: f32,
    samples: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PySurfacePoints> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let points = crate::execution::run(py, context, |context| {
        molframe::surface::surface_points(positions, radii, probe, samples, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(PySurfacePoints { points })
}

/// The accessible surface as points at `density` points per square ångström.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, density=4.0, context=None))]
fn surface_points_at_density(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    probe: f32,
    density: f32,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PySurfacePoints> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let points = crate::execution::run(py, context, |context| {
        molframe::surface::surface_points_at_density(positions, radii, probe, density, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(PySurfacePoints { points })
}

/// The solvent-excluded surface as a triangle mesh, on a grid of cell edge `resolution`.
///
/// `max_cells` and `max_workspace_bytes` raise the ceilings that keep the grid bounded.
#[pyfunction]
#[pyo3(signature = (
    coordinates,
    radii,
    *,
    probe=1.4,
    resolution=0.5,
    max_cells=None,
    max_workspace_bytes=None,
))]
fn solvent_excluded_surface(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    probe: f32,
    resolution: f32,
    max_cells: Option<usize>,
    max_workspace_bytes: Option<usize>,
) -> PyResult<crate::surface_mesh::PyMesh> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let options = grid(resolution, max_cells, max_workspace_bytes);
    let surface = py
        .detach(|| {
            molframe::surface::solvent_excluded_surface_with_options(
                positions, radii, probe, options,
            )
        })
        .map_err(crate::error::kernel)?;
    Ok(crate::surface_mesh::PyMesh::new(surface.indexed_mesh()))
}

/// The accessible areas that decide an interface, and the area it buries.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "BuriedSurface",
    frozen,
    skip_from_py_object,
    module = "molframe.surface"
)]
pub(crate) struct PyBuriedSurface {
    inner: molframe::surface::BuriedSurface,
}

#[pymethods]
impl PyBuriedSurface {
    /// Accessible area of the first group alone.
    #[getter]
    const fn first_alone(&self) -> f64 {
        self.inner.first_alone
    }

    /// Accessible area of the second group alone.
    #[getter]
    const fn second_alone(&self) -> f64 {
        self.inner.second_alone
    }

    /// Accessible area of both groups together.
    #[getter]
    const fn together(&self) -> f64 {
        self.inner.together
    }

    /// Area removed from the solvent by contact.
    #[getter]
    const fn buried(&self) -> f64 {
        self.inner.buried
    }

    fn __repr__(&self) -> String {
        format!("BuriedSurface(buried={:.3})", self.inner.buried)
    }
}

/// The surface buried between the atoms marked `in_first` and all the others.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, in_first, *, probe=1.4, points=960, context=None))]
fn buried_surface(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    in_first: &Bound<'_, PyArray1<bool>>,
    probe: f32,
    points: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyBuriedSurface> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let in_first = in_first.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let in_first = in_first.as_slice().map_err(|_| {
        crate::error::value("in_first must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    let inner = crate::execution::run(py, context, |context| {
        molframe::surface::buried_surface(positions, radii, probe, points, in_first, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(PyBuriedSurface { inner })
}

/// The solvent-excluded surface buried between two molecules.
///
/// `roles` names each atom's part: `first`, `second`, or `excluded` from every surface.
#[pyfunction]
#[pyo3(signature = (
    coordinates,
    radii,
    roles,
    *,
    probe=1.4,
    resolution=0.5,
    max_cells=None,
    max_workspace_bytes=None,
))]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
fn buried_solvent_excluded_surface(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    roles: Vec<String>,
    probe: f32,
    resolution: f32,
    max_cells: Option<usize>,
    max_workspace_bytes: Option<usize>,
) -> PyResult<PyBuriedSurface> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let roles = roles
        .iter()
        .map(|word| word.parse::<MoleculeRole>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(crate::error::kernel)?;
    let options = grid(resolution, max_cells, max_workspace_bytes);
    let inner = py
        .detach(|| {
            molframe::surface::buried_solvent_excluded_surface_with_options(
                positions, radii, &roles, probe, options,
            )
        })
        .map_err(crate::error::kernel)?;
    Ok(PyBuriedSurface { inner })
}

/// The distance from each atom to the nearest of `surface` points; infinity without any.
#[pyfunction]
#[pyo3(signature = (atoms, surface, *, cell_size=3.0))]
fn atom_depths<'py>(
    py: Python<'py>,
    atoms: &Bound<'py, PyArray2<f32>>,
    surface: &Bound<'py, PyArray2<f32>>,
    cell_size: f64,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let atoms = atoms.readonly();
    let surface = surface.readonly();
    let (atoms, surface) = (
        crate::bindings::coordinates(&atoms)?,
        crate::bindings::coordinates(&surface)?,
    );
    let depths = py
        .detach(|| molframe::surface::atom_depths(atoms, surface, AtomDepthOptions { cell_size }))
        .map_err(crate::error::kernel)?;
    Ok(depths.into_pyarray(py))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(sasa, module)?)?;
    module.add_function(wrap_pyfunction!(lee_richards, module)?)?;
    module.add_function(wrap_pyfunction!(cavities, module)?)?;
    module.add_function(wrap_pyfunction!(surface_points, module)?)?;
    module.add_function(wrap_pyfunction!(surface_points_at_density, module)?)?;
    module.add_function(wrap_pyfunction!(solvent_excluded_surface, module)?)?;
    module.add_function(wrap_pyfunction!(buried_surface, module)?)?;
    module.add_function(wrap_pyfunction!(buried_solvent_excluded_surface, module)?)?;
    module.add_function(wrap_pyfunction!(atom_depths, module)?)?;
    module.add_class::<PySurfacePoints>()?;
    module.add_class::<PyBuriedSurface>()?;
    crate::surface_mesh::register(module)
}
