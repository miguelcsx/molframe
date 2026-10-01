//! Mechanical adapters for reciprocal geometry and exact reflection symmetry.

use molframe::UnitCell;
use molframe::crystal::{CellTransform, ReflectionSymmetry, SymmetrySet, space_group_setting};
use pyo3::{exceptions::PyValueError, prelude::*};

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "UnitCell",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyUnitCell(CellTransform);

#[pymethods]
impl PyUnitCell {
    #[new]
    fn new(lengths: [f64; 3], angles: [f64; 3]) -> PyResult<Self> {
        CellTransform::new(&UnitCell { lengths, angles })
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn reciprocal_vector(&self, hkl: [i32; 3]) -> [f64; 3] {
        self.0.reciprocal_vector(hkl)
    }

    fn reciprocal_spacing_squared(&self, hkl: [i32; 3]) -> f64 {
        self.0.reciprocal_spacing_squared(hkl)
    }

    fn d_spacing(&self, hkl: [i32; 3]) -> f64 {
        self.0.d_spacing(hkl)
    }
}

#[derive(Clone, Debug)]
#[pyclass(
    name = "SpaceGroup",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PySpaceGroup(SymmetrySet);

#[pymethods]
impl PySpaceGroup {
    #[new]
    fn new(hall_number: u16) -> PyResult<Self> {
        space_group_setting(hall_number)
            .map(|setting| Self(setting.symmetry_set()))
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[getter]
    fn hall_symbol(&self) -> Option<&str> {
        self.0.hall.as_deref()
    }

    #[getter]
    fn hall_number(&self) -> Option<u16> {
        self.0.hall_number
    }

    fn reflection_symmetry(&self, hkl: [i32; 3]) -> PyResult<PyReflectionSymmetry> {
        self.0
            .reflection_symmetry(hkl)
            .map(PyReflectionSymmetry)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ReflectionSymmetry",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyReflectionSymmetry(ReflectionSymmetry);

#[pymethods]
impl PyReflectionSymmetry {
    #[getter]
    fn centric(&self) -> bool {
        self.0.centric
    }
    #[getter]
    fn systematically_absent(&self) -> bool {
        self.0.systematically_absent
    }
    #[getter]
    fn epsilon_factor(&self) -> usize {
        self.0.epsilon_factor
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyUnitCell>()?;
    module.add_class::<PySpaceGroup>()?;
    module.add_class::<PyReflectionSymmetry>()
}
