//! Owned reusable spatial state and safe Python adaptations for borrowed indexes.

use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use crate::spatial::{
    PyNeighborListOptions,
    spatial_arrays::{
        PyNeighborTable, SpatialBindingError, required_index_slice, validate_sorted_indices,
    },
};
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyType;

pyo3::create_exception!(pdbiox_spatial, SpatialError, PyValueError);

#[pyclass(name = "NeighborList", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyNeighborList(pdbiox::NeighborList);

#[pymethods]
impl PyNeighborList {
    #[classmethod]
    #[pyo3(signature = (positions, left, right, cutoff, options, *, cell=None))]
    fn build(
        _class: &Bound<'_, PyType>,
        py: Python<'_>,
        positions: PyReadonlyArray2<'_, f32>,
        left: PyReadonlyArray1<'_, u32>,
        right: PyReadonlyArray1<'_, u32>,
        cutoff: f32,
        options: PyNeighborListOptions,
        cell: Option<&PyUnitCell>,
    ) -> PyResult<Self> {
        let positions = borrowed_coordinates(&positions)?;
        let left = required_index_slice(&left)?;
        let right = required_index_slice(&right)?;
        let periodic = cell
            .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
            .transpose()
            .map_err(spatial_error)?;
        py.detach(|| {
            validate_sorted_indices(left, positions.len()).map_err(SpatialBindingError::Input)?;
            validate_sorted_indices(right, positions.len()).map_err(SpatialBindingError::Input)?;
            pdbiox::NeighborList::build_with_options(
                positions,
                left,
                right,
                cutoff,
                periodic.as_ref(),
                pdbiox::CoordinateGeneration::INITIAL,
                options.0,
            )
            .map_err(SpatialBindingError::Spatial)
        })
        .map(Self)
        .map_err(SpatialBindingError::into_pyerr)
    }

    fn generation(&self) -> u64 {
        self.0.generation().get()
    }

    fn is_current(&self, generation: u64) -> bool {
        self.0.generation().get() == generation
    }

    #[pyo3(signature = (positions, *, cell=None))]
    fn can_reuse(
        &self,
        py: Python<'_>,
        positions: PyReadonlyArray2<'_, f32>,
        cell: Option<&PyUnitCell>,
    ) -> PyResult<bool> {
        let positions = borrowed_coordinates(&positions)?;
        let periodic = cell
            .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
            .transpose()
            .map_err(spatial_error)?;
        Ok(py.detach(|| self.0.can_reuse(positions, periodic.as_ref())))
    }

    #[pyo3(signature = (positions, cutoff, *, cell=None))]
    fn pairs(
        &self,
        py: Python<'_>,
        positions: PyReadonlyArray2<'_, f32>,
        cutoff: f32,
        cell: Option<&PyUnitCell>,
    ) -> PyResult<PyNeighborTable> {
        let positions = borrowed_coordinates(&positions)?;
        let periodic = cell
            .map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
            .transpose()
            .map_err(spatial_error)?;
        py.detach(|| {
            let pairs = self
                .0
                .pairs(positions, cutoff, periodic.as_ref())
                .map_err(SpatialBindingError::Spatial)?;
            PyNeighborTable::from_pairs(pairs).map_err(SpatialBindingError::Allocation)
        })
        .map_err(SpatialBindingError::into_pyerr)
    }
}

pub(crate) fn spatial_error(error: impl std::fmt::Display) -> PyErr {
    SpatialError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("SpatialError", module.py().get_type::<SpatialError>())?;
    module.add_class::<PyNeighborList>()?;
    Ok(())
}
