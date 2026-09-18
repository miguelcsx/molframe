//! Vectorized native neighbour queries over explicit coordinate arrays.

use crate::graph::PySpatialBackend;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[path = "spatial_arrays.rs"]
pub(crate) mod spatial_arrays;
pub(crate) use spatial_arrays::PyNeighborTable;

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

#[pymethods]
impl PyNeighborPair {
    #[new]
    fn new(first: u32, second: u32, distance: f32) -> PyResult<Self> {
        if first >= second {
            return Err(PyValueError::new_err(
                "neighbor pairs require distinct ascending atom indices",
            ));
        }
        if !distance.is_finite() || distance < 0.0 {
            return Err(PyValueError::new_err(
                "neighbor-pair distance must be finite and non-negative",
            ));
        }
        Ok(Self {
            first,
            second,
            distance,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "NeighborPair(first={}, second={}, distance={})",
            self.first, self.second, self.distance
        )
    }
}

#[pyclass(name = "AutoBackendProfile", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAutoBackendProfile(molframe::AutoBackendProfile);

#[pymethods]
impl PyAutoBackendProfile {
    #[new]
    fn new(brute_force_pair_limit: usize, kd_target_minimum: usize, kd_query_ratio: usize) -> Self {
        Self(molframe::AutoBackendProfile {
            brute_force_pair_limit,
            kd_target_minimum,
            kd_query_ratio,
        })
    }

    #[staticmethod]
    fn balanced() -> Self {
        Self(molframe::AutoBackendProfile::BALANCED)
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
        Self::from(molframe::NeighborSkinProfile::BALANCED)
    }
}

impl PyNeighborSkinProfile {
    const fn inner(self) -> molframe::NeighborSkinProfile {
        molframe::NeighborSkinProfile {
            cutoff_ratio: self.cutoff_ratio,
            minimum: self.minimum,
        }
    }
}

impl From<molframe::NeighborSkinProfile> for PyNeighborSkinProfile {
    fn from(value: molframe::NeighborSkinProfile) -> Self {
        Self {
            cutoff_ratio: value.cutoff_ratio,
            minimum: value.minimum,
        }
    }
}

#[pyclass(name = "CellGridOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCellGridOptions(molframe::CellGridOptions);

#[pymethods]
impl PyCellGridOptions {
    #[new]
    fn new(maximum_cell_count: usize, edge_growth_factor: f32) -> Self {
        Self(molframe::CellGridOptions {
            maximum_cell_count,
            edge_growth_factor,
        })
    }

    #[staticmethod]
    fn memory_balanced() -> Self {
        Self(molframe::CellGridOptions::MEMORY_BALANCED)
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

#[pyclass(name = "NeighborListOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNeighborListOptions(pub(crate) molframe::NeighborListOptions);

#[pymethods]
impl PyNeighborListOptions {
    #[new]
    #[pyo3(signature = (skin, *, cell_grid=None))]
    fn new(skin: f32, cell_grid: Option<PyCellGridOptions>) -> Self {
        let mut inner = molframe::NeighborListOptions::with_skin(skin);
        if let Some(cell_grid) = cell_grid {
            inner.cell_grid = cell_grid.0;
        }
        Self(inner)
    }

    #[staticmethod]
    fn with_skin(skin: f32) -> Self {
        Self(molframe::NeighborListOptions::with_skin(skin))
    }

    #[getter]
    fn skin(&self) -> f32 {
        self.0.skin
    }

    #[getter]
    fn cell_grid(&self) -> PyCellGridOptions {
        PyCellGridOptions(self.0.cell_grid)
    }

    fn validate(&self) -> PyResult<()> {
        self.0.validate().map_err(value_error)
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
            maximum_image_count: molframe::KdPeriodicOptions::BALANCED.maximum_image_count,
        }
    }
}

impl PyKdPeriodicOptions {
    pub(crate) const fn inner(self) -> molframe::KdPeriodicOptions {
        molframe::KdPeriodicOptions {
            maximum_image_count: self.maximum_image_count,
        }
    }
}

#[pyclass(name = "SpatialSearchOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySpatialSearchOptions(molframe::SpatialSearchOptions);

impl PySpatialSearchOptions {
    pub(crate) const fn from_native(value: molframe::SpatialSearchOptions) -> Self {
        Self(value)
    }

    pub(crate) const fn balanced_value() -> Self {
        Self::from_native(molframe::SpatialSearchOptions::BALANCED)
    }

    pub(crate) const fn inner(self) -> molframe::SpatialSearchOptions {
        self.0
    }
}

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
        let balanced = molframe::SpatialSearchOptions::BALANCED;
        Self(molframe::SpatialSearchOptions {
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
        Self(molframe::SpatialSearchOptions::BALANCED)
    }

    fn plan(&self, left_count: usize, right_count: usize, cutoff: f32) -> PyResult<PySpatialPlan> {
        self.0
            .plan(left_count, right_count, cutoff)
            .map(PySpatialPlan)
            .map_err(value_error)
    }
}

#[pyclass(name = "SpatialPlan", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySpatialPlan(molframe::SpatialPlan);

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

pub(crate) fn normalize_selection(
    indices: Option<Vec<u32>>,
    count: usize,
) -> Result<molframe::AtomSelection, String> {
    let Some(mut indices) = indices else {
        return u32::try_from(count)
            .map(molframe::AtomSelection::All)
            .map_err(|_| "coordinate count exceeds u32 capacity".to_owned());
    };
    indices.sort_unstable();
    indices.dedup();
    if indices.last().is_some_and(|index| *index as usize >= count) {
        return Err("atom index is outside the coordinate array".to_owned());
    }
    Ok(molframe::AtomSelection::from_sorted(indices))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::spatial_index::register(module)?;
    module.add_class::<crate::spatial_periodic::PyPeriodicImage>()?;
    module.add_class::<crate::spatial_periodic::PyPeriodicBox>()?;
    module.add_class::<crate::spatial_periodic::PySpatialOption>()?;
    module.add_class::<PyNeighborPair>()?;
    module.add_class::<PyAutoBackendProfile>()?;
    module.add_class::<PyNeighborSkinProfile>()?;
    module.add_class::<PyCellGridOptions>()?;
    module.add_class::<PyNeighborListOptions>()?;
    module.add_class::<PyKdPeriodicOptions>()?;
    module.add_class::<PySpatialSearchOptions>()?;
    module.add_class::<PySpatialPlan>()?;
    spatial_arrays::register(module)?;
    super::count::register(module)?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    crate::spatial_index::spatial_error(error)
}

#[cfg(test)]
#[path = "spatial_tests.rs"]
mod tests;
