//! Governed physical fields of one structure: weighted densities, pore
//! profiles and exposed-surface contacts.
//!
//! Every binning and every threshold is a required keyword; weights are a
//! name or an explicit per-atom array.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use crate::governed::run;
use crate::governed_physical::{axis_of, weights_of};
use crate::policy::{PyAnalysisPolicy, policy_of};
use crate::table::TableBuilder;
use molframe::analysis::{
    DensityGridSpec, LinearDensityOptions, PoreProfileOptions, SurfaceContactOptions,
};
use numpy::{PyArray3, PyArrayMethods, ToPyArray};
use pyo3::prelude::*;

/// A weighted density profile along one Cartesian axis.
///
/// `weights` is `"count"`, `"mass"` or one value per atom (charge, occupancy,
/// …): no unit is inferred. Bins are half-open `[lower, upper)`. Columns:
/// `lower`, `upper`, `weight`, `density` (weight per ångström).
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    axis,
    minimum,
    maximum,
    bins,
    weights=None,
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn linear_density(
    py: Python<'_>,
    structure: &PyStructure,
    axis: &str,
    minimum: f32,
    maximum: f32,
    bins: usize,
    weights: Option<&Bound<'_, PyAny>>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let weights = weights_of(structure, weights)?;
    let kernel = molframe::analysis::linear_density_kernel(
        &weights,
        LinearDensityOptions {
            axis: axis_of(axis)?,
            minimum,
            maximum,
            bins,
        },
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let lower: Vec<f32> = rows.iter().map(|row| row.lower).collect();
            let upper: Vec<f32> = rows.iter().map(|row| row.upper).collect();
            let weight: Vec<f64> = rows.iter().map(|row| row.weight).collect();
            let density: Vec<f64> = rows.iter().map(|row| row.density).collect();
            let table = TableBuilder::new(py, rows.len())
                .single("lower", &lower)
                .single("upper", &upper)
                .double("weight", &weight)
                .double("density", &density)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// A weighted density on a rectilinear grid.
#[pyclass(
    name = "DensityGrid",
    module = "molframe.analysis",
    frozen,
    skip_from_py_object
)]
pub(crate) struct PyDensityGrid {
    #[pyo3(get)]
    origin: [f32; 3],
    #[pyo3(get)]
    spacing: [f32; 3],
    #[pyo3(get)]
    shape: [usize; 3],
    #[pyo3(get)]
    excluded_weight: f64,
    values: Vec<f64>,
}

#[pymethods]
impl PyDensityGrid {
    /// Weight per unit volume, shape `(nx, ny, nz)`; voxel `(i, j, k)` spans
    /// `origin + spacing * (i, j, k)` to the next voxel corner.
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let array =
            self.values
                .to_pyarray(py)
                .reshape([self.shape[0], self.shape[1], self.shape[2]])?;
        array.readwrite().make_nonwriteable();
        Ok(array)
    }

    fn __repr__(&self) -> String {
        format!(
            "DensityGrid(shape={:?}, excluded_weight={})",
            self.shape, self.excluded_weight
        )
    }
}

/// A weighted density on the grid with lower corner `origin`, voxel edges
/// `spacing` (Å) and `shape` voxels; weight outside the grid is reported as
/// `excluded_weight`, never dropped silently.
#[pyfunction]
#[pyo3(signature = (structure, *, origin, spacing, shape, weights=None, policy=None, context=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn density_map(
    py: Python<'_>,
    structure: &PyStructure,
    origin: [f32; 3],
    spacing: [f32; 3],
    shape: [usize; 3],
    weights: Option<&Bound<'_, PyAny>>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let weights = weights_of(structure, weights)?;
    let kernel = molframe::analysis::density_map_kernel(
        &weights,
        DensityGridSpec {
            origin,
            spacing,
            shape,
        },
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, grid| {
            Ok(Py::new(
                py,
                PyDensityGrid {
                    origin: grid.spec.origin,
                    spacing: grid.spec.spacing,
                    shape: grid.spec.shape,
                    excluded_weight: grid.excluded_weight,
                    values: grid.density,
                },
            )?
            .into_any())
        },
        context,
    )
}

/// The largest probe radius that fits at each of `samples` slices along an axis.
///
/// The axis passes through `axis_origin` along `axis_direction`; slices are
/// evenly spaced from `start` to `end` along it. At each slice a transverse
/// grid of spacing `grid_spacing` out to `search_radius` is searched for the
/// centre with the greatest clearance, from which `probe_radius` is
/// subtracted. Columns: `axial_coordinate`, `centre_x`, `centre_y`,
/// `centre_z`, `radius`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    axis_origin,
    axis_direction,
    start,
    end,
    samples,
    search_radius,
    grid_spacing,
    probe_radius=0.0,
    radii="bondi",
    memory_limit_bytes=PoreProfileOptions::DEFAULT_MEMORY_LIMIT_BYTES,
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn pore_profile(
    py: Python<'_>,
    structure: &PyStructure,
    axis_origin: [f32; 3],
    axis_direction: [f32; 3],
    start: f32,
    end: f32,
    samples: usize,
    search_radius: f32,
    grid_spacing: f32,
    probe_radius: f32,
    radii: &str,
    memory_limit_bytes: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let radii = molframe::chemistry::atom_radii(
        structure.inner.engine(),
        radii.parse().map_err(crate::error::kernel)?,
    );
    let kernel = molframe::analysis::pore_profile_kernel(
        &radii,
        PoreProfileOptions::new(
            (axis_origin, axis_direction),
            start,
            end,
            samples,
            search_radius,
            grid_spacing,
            probe_radius,
        )
        .with_memory_limit(memory_limit_bytes),
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let axial: Vec<f32> = rows.iter().map(|row| row.axial_coordinate).collect();
            let axis =
                |index: usize| -> Vec<f32> { rows.iter().map(|row| row.centre[index]).collect() };
            let radius: Vec<f32> = rows.iter().map(|row| row.radius).collect();
            let table = TableBuilder::new(py, rows.len())
                .single("axial_coordinate", &axial)
                .single("centre_x", &axis(0))
                .single("centre_y", &axis(1))
                .single("centre_z", &axis(2))
                .single("radius", &radius)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Atom pairs within `tolerance` ångström of touching whose facing surfaces
/// are both exposed.
///
/// A pair qualifies when the distance is at most the sum of the two radii plus
/// `tolerance` and each atom has at least `minimum_area` Å² of solvent-accessible
/// surface (probe `probe`, `surface_density` points per Å²). Columns: `first`,
/// `second`, `distance`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    tolerance,
    probe,
    surface_density,
    minimum_area,
    radii="bondi",
    backend="auto",
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn surface_contacts(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    probe: f32,
    surface_density: f32,
    minimum_area: f32,
    radii: &str,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let radii = molframe::chemistry::atom_radii(
        structure.inner.engine(),
        radii.parse().map_err(crate::error::kernel)?,
    );
    let kernel = molframe::analysis::surface_contacts_kernel(
        &radii,
        SurfaceContactOptions {
            tolerance,
            probe,
            surface_density,
            minimum_area,
            backend: crate::backend::parse(backend)?,
        },
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let first: Vec<u32> = rows.iter().map(|row| row.first.get()).collect();
            let second: Vec<u32> = rows.iter().map(|row| row.second.get()).collect();
            let distance: Vec<f32> = rows.iter().map(|row| row.distance).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("first", &first)
                .indices("second", &second)
                .single("distance", &distance)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Atom contacts under the policy's `contact_def` and `vdw_radii`.
///
/// With `distance:<tolerance>` two atoms are in contact when their distance is at
/// most the sum of their radii plus the tolerance; an atom whose element has no
/// radius in the set is missing coverage. With `surface:<probe>` they are in contact
/// when their expanded surfaces touch and both keep an exposed patch; the three
/// `surface_*` keywords are how that test samples the surface and apply only to it.
/// Columns: `first`, `second`, `distance`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    backend="auto",
    surface_tolerance=0.2,
    surface_density=4.0,
    surface_minimum_area=0.25,
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn contacts_by_definition(
    py: Python<'_>,
    structure: &PyStructure,
    backend: &str,
    surface_tolerance: f32,
    surface_density: f32,
    surface_minimum_area: f32,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::analysis::definition_contacts_kernel(
        crate::backend::parse(backend)?,
        molframe::analysis::SurfaceSampling {
            tolerance: surface_tolerance,
            density: surface_density,
            minimum_area: surface_minimum_area,
        },
    );
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
            let first: Vec<u32> = table.first().iter().map(|atom| atom.get()).collect();
            let second: Vec<u32> = table.second().iter().map(|atom| atom.get()).collect();
            let table = TableBuilder::new(py, table.len())
                .indices("first", &first)
                .indices("second", &second)
                .single("distance", table.distances())
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDensityGrid>()?;
    module.add_function(wrap_pyfunction!(linear_density, module)?)?;
    module.add_function(wrap_pyfunction!(density_map, module)?)?;
    module.add_function(wrap_pyfunction!(pore_profile, module)?)?;
    module.add_function(wrap_pyfunction!(surface_contacts, module)?)?;
    module.add_function(wrap_pyfunction!(contacts_by_definition, module)?)
}
