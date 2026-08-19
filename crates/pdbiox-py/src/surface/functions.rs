//! Python entry points for the reusable surface kernels.

use crate::geometry::borrowed_coordinates;
use crate::intrinsic::PySurfaceGridOptions;
use crate::surface_types::{
    PyAtomContactArea, PyAtomDepthOptions, PyBuriedSurface, PyCavity, PyExcludedSurfacePoint,
    PyIndexedSurfaceMesh, PyMoleculeRole, PySolventExcludedSurface, PySurfaceCurvature,
    PySurfaceDistances, PySurfacePoint,
};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn fibonacci_sphere(count: u16) -> Vec<[f64; 3]> {
    pdbiox::surface::fibonacci_sphere(count)
}

#[pyfunction]
pub(crate) fn shrake_rupley(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    points: u16,
) -> PyResult<Vec<f64>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::shrake_rupley(positions, radii, probe, points))
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn lee_richards(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    slices: u16,
) -> PyResult<Vec<f64>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::lee_richards(positions, radii, probe, slices))
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
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::shrake_rupley(positions, radii, probe_radius, sample_points))
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
) -> PyResult<PyBuriedSurface> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    let first = first.as_slice()?;
    py.detach(move || {
        pdbiox::surface::buried_surface(positions, radii, probe_radius, sample_points, first)
    })
    .map(Into::into)
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
    let atoms = borrowed_coordinates(&atoms)?;
    let surface = borrowed_coordinates(&surface)?;
    py.detach(move || pdbiox::surface::atom_depths(atoms, surface, options.0))
        .map(|values| values.into_pyarray(py))
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn cavities(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<Vec<PyCavity>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::cavities_with_options(positions, radii, probe, options.0))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn surface_points(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    samples: u16,
) -> PyResult<Vec<PySurfacePoint>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::surface_points(positions, radii, probe, samples))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn surface_points_at_density(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    density: f32,
) -> PyResult<Vec<PySurfacePoint>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::surface_points_at_density(positions, radii, probe, density))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn atom_contact_areas(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    density: f32,
) -> PyResult<Vec<PyAtomContactArea>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || pdbiox::surface::atom_contact_areas(positions, radii, probe, density))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn surface_points_excluding_pairs(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    density: f32,
    pairs: Vec<[usize; 2]>,
) -> PyResult<Vec<PyExcludedSurfacePoint>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    let pairs = pairs
        .into_iter()
        .map(|pair| (pair[0], pair[1]))
        .collect::<Vec<_>>();
    py.detach(move || {
        pdbiox::surface::surface_points_excluding_pairs(positions, radii, probe, density, &pairs)
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn buried_solvent_excluded_surface(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    roles: Vec<PyMoleculeRole>,
    probe: f32,
    resolution: f32,
) -> PyResult<PyBuriedSurface> {
    buried_solvent_excluded_surface_with_options_impl(
        py,
        positions,
        radii,
        roles,
        probe,
        PySurfaceGridOptions::standard(resolution),
    )
}

#[pyfunction]
pub(crate) fn buried_solvent_excluded_surface_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    roles: Vec<PyMoleculeRole>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<PyBuriedSurface> {
    buried_solvent_excluded_surface_with_options_impl(py, positions, radii, roles, probe, options)
}

fn buried_solvent_excluded_surface_with_options_impl(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    roles: Vec<PyMoleculeRole>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<PyBuriedSurface> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    let roles = roles.into_iter().map(Into::into).collect::<Vec<_>>();
    py.detach(move || {
        pdbiox::surface::buried_solvent_excluded_surface_with_options(
            positions, radii, &roles, probe, options.0,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn solvent_excluded_surface(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    resolution: f32,
) -> PyResult<PySolventExcludedSurface> {
    solvent_excluded_surface_with_options_impl(
        py,
        positions,
        radii,
        probe,
        PySurfaceGridOptions::standard(resolution),
    )
}

#[pyfunction]
pub(crate) fn solvent_excluded_surface_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<PySolventExcludedSurface> {
    solvent_excluded_surface_with_options_impl(py, positions, radii, probe, options)
}

fn solvent_excluded_surface_with_options_impl(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    options: PySurfaceGridOptions,
) -> PyResult<PySolventExcludedSurface> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    py.detach(move || {
        pdbiox::surface::solvent_excluded_surface_with_options(positions, radii, probe, options.0)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn edge_geodesic_distances(
    mesh: &PyIndexedSurfaceMesh,
    source: u32,
) -> PyResult<PySurfaceDistances> {
    pdbiox::surface::edge_geodesic_distances(&mesh.native, source)
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn surface_patch(
    mesh: &PyIndexedSurfaceMesh,
    source: u32,
    radius: f64,
) -> PyResult<Vec<u32>> {
    pdbiox::surface::surface_patch(&mesh.native, source, radius)
        .map(<[u32]>::into_vec)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn surface_curvatures(mesh: &PyIndexedSurfaceMesh) -> PyResult<Vec<PySurfaceCurvature>> {
    pdbiox::surface::surface_curvatures(&mesh.native)
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn write_obj(path: &str, mesh: &PyIndexedSurfaceMesh) -> PyResult<()> {
    pdbiox::surface::write_obj(path, &mesh.native).map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(fibonacci_sphere, module)?)?;
    module.add_function(wrap_pyfunction!(shrake_rupley, module)?)?;
    module.add_function(wrap_pyfunction!(lee_richards, module)?)?;
    module.add_function(wrap_pyfunction!(solvent_accessible_surface, module)?)?;
    module.add_function(wrap_pyfunction!(buried_surface, module)?)?;
    module.add_function(wrap_pyfunction!(atom_depths, module)?)?;
    module.add_function(wrap_pyfunction!(cavities, module)?)?;
    module.add_function(wrap_pyfunction!(surface_points, module)?)?;
    module.add_function(wrap_pyfunction!(surface_points_at_density, module)?)?;
    module.add_function(wrap_pyfunction!(atom_contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(surface_points_excluding_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(buried_solvent_excluded_surface, module)?)?;
    module.add_function(wrap_pyfunction!(
        buried_solvent_excluded_surface_with_options,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(solvent_excluded_surface, module)?)?;
    module.add_function(wrap_pyfunction!(
        solvent_excluded_surface_with_options,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(edge_geodesic_distances, module)?)?;
    module.add_function(wrap_pyfunction!(surface_patch, module)?)?;
    module.add_function(wrap_pyfunction!(surface_curvatures, module)?)?;
    module.add_function(wrap_pyfunction!(write_obj, module)?)?;
    Ok(())
}
