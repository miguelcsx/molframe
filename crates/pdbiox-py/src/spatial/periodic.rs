//! Direct periodic-geometry projections over the native spatial kernel.

use crate::crystallography::PyUnitCell;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "PeriodicImage", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeriodicImage {
    #[pyo3(get)]
    displacement: [f32; 3],
    #[pyo3(get)]
    lattice_shift: [i64; 3],
}

impl From<pdbiox::PeriodicImage> for PyPeriodicImage {
    fn from(value: pdbiox::PeriodicImage) -> Self {
        Self {
            displacement: value.displacement,
            lattice_shift: value.lattice_shift,
        }
    }
}

#[pyclass(name = "PeriodicBox", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeriodicBox(pdbiox::PeriodicBox);

#[pymethods]
impl PyPeriodicBox {
    #[new]
    fn new(cell: &PyUnitCell) -> PyResult<Self> {
        pdbiox::PeriodicBox::from_cell(cell.cell)
            .map(Self)
            .map_err(value_error)
    }

    fn displacement(&self, left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
        self.0.displacement(left, right)
    }

    fn minimum_image(&self, left: [f32; 3], right: [f32; 3]) -> PyPeriodicImage {
        self.0.minimum_image(left, right).into()
    }

    fn distance_squared(&self, left: [f32; 3], right: [f32; 3]) -> f32 {
        self.0.distance_squared(left, right)
    }

    fn fractional(&self, position: [f32; 3]) -> [f64; 3] {
        self.0.fractional(position)
    }

    fn cartesian(&self, fractional: [f64; 3]) -> [f32; 3] {
        self.0.cartesian(fractional)
    }

    fn wrap(&self, position: [f32; 3]) -> [f32; 3] {
        self.0.wrap(position)
    }

    fn interpolate(&self, left: [f32; 3], right: [f32; 3], amount: f64) -> Option<[f32; 3]> {
        self.0.interpolate(left, right, amount)
    }
}

#[pyclass(name = "SpatialOption", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySpatialOption {
    BruteForcePairLimit,
    KdTargetMinimum,
    KdQueryRatio,
    PeriodicBackend,
    NeighborSkinRatio,
    NeighborSkinMinimum,
    MaximumCellCount,
    CellGrowthFactor,
    KdPeriodicImageLimit,
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
