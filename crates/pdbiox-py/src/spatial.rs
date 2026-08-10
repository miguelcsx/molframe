//! Vectorized native neighbour queries over explicit coordinate arrays.

use crate::crystallography::PyUnitCell;
use crate::geometry::coordinates;
use crate::graph::PySpatialBackend;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "NeighborPair", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNeighborPair {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
    #[pyo3(get)]
    distance: f32,
}

#[pyclass(name = "AutoBackendProfile", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAutoBackendProfile(pdbiox::AutoBackendProfile);

#[pymethods]
impl PyAutoBackendProfile {
    #[new]
    fn new(
        brute_force_pair_limit: usize,
        kd_target_minimum: usize,
        kd_query_ratio: usize,
        periodic_backend: PySpatialBackend,
    ) -> Self {
        Self(pdbiox::AutoBackendProfile {
            brute_force_pair_limit,
            kd_target_minimum,
            kd_query_ratio,
            periodic_backend: periodic_backend.into(),
        })
    }

    #[staticmethod]
    fn balanced() -> Self {
        Self(pdbiox::AutoBackendProfile::BALANCED)
    }

    #[getter]
    fn brute_force_pair_limit(&self) -> usize {
        self.0.brute_force_pair_limit
    }
    #[getter]
    fn kd_target_minimum(&self) -> usize {
        self.0.kd_target_minimum
    }
    #[getter]
    fn kd_query_ratio(&self) -> usize {
        self.0.kd_query_ratio
    }
    #[getter]
    fn periodic_backend(&self) -> PySpatialBackend {
        self.0.periodic_backend.into()
    }
}

#[pyclass(name = "NeighborSkinProfile", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNeighborSkinProfile {
    #[pyo3(get)]
    cutoff_ratio: f32,
    #[pyo3(get)]
    minimum: f32,
}

#[pymethods]
impl PyNeighborSkinProfile {
    #[new]
    fn new(cutoff_ratio: f32, minimum: f32) -> Self {
        Self {
            cutoff_ratio,
            minimum,
        }
    }

    #[staticmethod]
    fn balanced() -> Self {
        Self::from(pdbiox::NeighborSkinProfile::BALANCED)
    }
}

impl PyNeighborSkinProfile {
    const fn inner(self) -> pdbiox::NeighborSkinProfile {
        pdbiox::NeighborSkinProfile {
            cutoff_ratio: self.cutoff_ratio,
            minimum: self.minimum,
        }
    }
}

impl From<pdbiox::NeighborSkinProfile> for PyNeighborSkinProfile {
    fn from(value: pdbiox::NeighborSkinProfile) -> Self {
        Self {
            cutoff_ratio: value.cutoff_ratio,
            minimum: value.minimum,
        }
    }
}

#[pyclass(name = "CellGridOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCellGridOptions(pdbiox::CellGridOptions);

#[pymethods]
impl PyCellGridOptions {
    #[new]
    fn new(maximum_cell_count: usize, edge_growth_factor: f32) -> Self {
        Self(pdbiox::CellGridOptions {
            maximum_cell_count,
            edge_growth_factor,
        })
    }

    #[staticmethod]
    fn memory_balanced() -> Self {
        Self(pdbiox::CellGridOptions::MEMORY_BALANCED)
    }

    #[getter]
    fn maximum_cell_count(&self) -> usize {
        self.0.maximum_cell_count
    }
    #[getter]
    fn edge_growth_factor(&self) -> f32 {
        self.0.edge_growth_factor
    }
}

#[pyclass(name = "KdPeriodicOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyKdPeriodicOptions {
    #[pyo3(get)]
    maximum_image_count: usize,
}

#[pymethods]
impl PyKdPeriodicOptions {
    #[new]
    fn new(maximum_image_count: usize) -> Self {
        Self {
            maximum_image_count,
        }
    }

    #[staticmethod]
    fn balanced() -> Self {
        Self {
            maximum_image_count: pdbiox::KdPeriodicOptions::BALANCED.maximum_image_count,
        }
    }
}

impl PyKdPeriodicOptions {
    const fn inner(self) -> pdbiox::KdPeriodicOptions {
        pdbiox::KdPeriodicOptions {
            maximum_image_count: self.maximum_image_count,
        }
    }
}

#[pyclass(name = "SpatialSearchOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySpatialSearchOptions(pdbiox::SpatialSearchOptions);

#[pymethods]
impl PySpatialSearchOptions {
    #[new]
    #[pyo3(signature = (backend=PySpatialBackend::Auto, *, automatic=None, neighbor_skin=None, cell_grid=None, kd_periodic=None))]
    fn new(
        backend: PySpatialBackend,
        automatic: Option<PyAutoBackendProfile>,
        neighbor_skin: Option<PyNeighborSkinProfile>,
        cell_grid: Option<PyCellGridOptions>,
        kd_periodic: Option<PyKdPeriodicOptions>,
    ) -> Self {
        let balanced = pdbiox::SpatialSearchOptions::BALANCED;
        Self(pdbiox::SpatialSearchOptions {
            backend: backend.into(),
            automatic: automatic.map_or(balanced.automatic, |value| value.0),
            neighbor_skin: neighbor_skin
                .map_or(balanced.neighbor_skin, PyNeighborSkinProfile::inner),
            cell_grid: cell_grid.map_or(balanced.cell_grid, |value| value.0),
            kd_periodic: kd_periodic.map_or(balanced.kd_periodic, PyKdPeriodicOptions::inner),
        })
    }

    #[staticmethod]
    fn balanced() -> Self {
        Self(pdbiox::SpatialSearchOptions::BALANCED)
    }

    fn plan(
        &self,
        left_count: usize,
        right_count: usize,
        periodic: bool,
        cutoff: f32,
    ) -> PyResult<PySpatialPlan> {
        self.0
            .plan(left_count, right_count, periodic, cutoff)
            .map(PySpatialPlan)
            .map_err(value_error)
    }
}

#[pyclass(name = "SpatialPlan", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySpatialPlan(pdbiox::SpatialPlan);

#[pymethods]
impl PySpatialPlan {
    #[getter]
    fn requested_backend(&self) -> PySpatialBackend {
        self.0.requested_backend.into()
    }
    #[getter]
    fn backend(&self) -> PySpatialBackend {
        self.0.backend.into()
    }
    #[getter]
    fn neighbor_skin(&self) -> Option<f32> {
        self.0.neighbor_skin
    }
    #[getter]
    fn options(&self) -> PySpatialSearchOptions {
        PySpatialSearchOptions(self.0.options)
    }
}

#[pyfunction]
#[pyo3(signature = (positions, cutoff, *, left=None, right=None, backend=PySpatialBackend::Auto, cell=None))]
pub(crate) fn neighbor_pairs(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    cutoff: f32,
    left: Option<Vec<u32>>,
    right: Option<Vec<u32>>,
    backend: PySpatialBackend,
    cell: Option<&PyUnitCell>,
) -> PyResult<Vec<PyNeighborPair>> {
    let positions = coordinates(positions)?;
    let left = selection(left, positions.len())?;
    let right = selection(right, positions.len())?;
    let periodic = cell
        .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
        .transpose()
        .map_err(value_error)?;
    py.detach(move || {
        pdbiox::pairs_within(
            &positions,
            &left,
            &right,
            cutoff,
            backend.into(),
            periodic.as_ref(),
        )
    })
    .map(|values| values.into_iter().map(PyNeighborPair::from).collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (positions, targets, cutoff, *, query=None, backend=PySpatialBackend::Auto, cell=None))]
pub(crate) fn atoms_within(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    targets: Vec<u32>,
    cutoff: f32,
    query: Option<Vec<u32>>,
    backend: PySpatialBackend,
    cell: Option<&PyUnitCell>,
) -> PyResult<Vec<u32>> {
    let positions = coordinates(positions)?;
    let query = selection(query, positions.len())?;
    let targets = selection(Some(targets), positions.len())?;
    let periodic = cell
        .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
        .transpose()
        .map_err(value_error)?;
    py.detach(move || {
        pdbiox::within(
            &positions,
            &query,
            &targets,
            cutoff,
            backend.into(),
            periodic.as_ref(),
        )
    })
    .map(|selection| selection.into_iter().collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (positions, cutoff, options, *, left=None, right=None, cell=None))]
pub(crate) fn neighbor_pairs_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    cutoff: f32,
    options: PySpatialSearchOptions,
    left: Option<Vec<u32>>,
    right: Option<Vec<u32>>,
    cell: Option<&PyUnitCell>,
) -> PyResult<Vec<PyNeighborPair>> {
    let positions = coordinates(positions)?;
    let left = selection(left, positions.len())?;
    let right = selection(right, positions.len())?;
    let periodic = cell
        .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
        .transpose()
        .map_err(value_error)?;
    py.detach(move || {
        pdbiox::pairs_within_with_options(
            &positions,
            &left,
            &right,
            cutoff,
            options.0,
            periodic.as_ref(),
        )
    })
    .map(|values| values.into_iter().map(PyNeighborPair::from).collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (positions, targets, cutoff, options, *, query=None, cell=None))]
pub(crate) fn atoms_within_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    targets: Vec<u32>,
    cutoff: f32,
    options: PySpatialSearchOptions,
    query: Option<Vec<u32>>,
    cell: Option<&PyUnitCell>,
) -> PyResult<Vec<u32>> {
    let positions = coordinates(positions)?;
    let query = selection(query, positions.len())?;
    let targets = selection(Some(targets), positions.len())?;
    let periodic = cell
        .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
        .transpose()
        .map_err(value_error)?;
    py.detach(move || {
        pdbiox::within_with_options(
            &positions,
            &query,
            &targets,
            cutoff,
            options.0,
            periodic.as_ref(),
        )
    })
    .map(|selection| selection.into_iter().collect())
    .map_err(value_error)
}

fn selection(indices: Option<Vec<u32>>, count: usize) -> PyResult<pdbiox::AtomSelection> {
    let Some(mut indices) = indices else {
        return u32::try_from(count)
            .map(pdbiox::AtomSelection::All)
            .map_err(|_| PyValueError::new_err("coordinate count exceeds u32 capacity"));
    };
    indices.sort_unstable();
    indices.dedup();
    if indices.last().is_some_and(|index| *index as usize >= count) {
        return Err(PyValueError::new_err(
            "atom index is outside the coordinate array",
        ));
    }
    Ok(pdbiox::AtomSelection::from_sorted(indices))
}

impl From<pdbiox::NeighborPair> for PyNeighborPair {
    fn from(value: pdbiox::NeighborPair) -> Self {
        Self {
            first: value.first,
            second: value.second,
            distance: value.distance_squared.sqrt(),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNeighborPair>()?;
    module.add_class::<PyAutoBackendProfile>()?;
    module.add_class::<PyNeighborSkinProfile>()?;
    module.add_class::<PyCellGridOptions>()?;
    module.add_class::<PyKdPeriodicOptions>()?;
    module.add_class::<PySpatialSearchOptions>()?;
    module.add_class::<PySpatialPlan>()?;
    module.add_function(wrap_pyfunction!(neighbor_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(atoms_within, module)?)?;
    module.add_function(wrap_pyfunction!(neighbor_pairs_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(atoms_within_with_options, module)?)?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
