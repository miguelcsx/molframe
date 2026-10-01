//! Mechanical adapters for reciprocal geometry and exact reflection symmetry.

use crate::bindings::PyStructure;
use molframe::UnitCell;
use molframe::crystal::{CellTransform, ReflectionSymmetry, SymmetrySet, space_group_setting};
use numpy::{Complex64, PyArray1, PyArray2, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
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

/// X-ray structure factors of the first model, one complex value per `(h, k, l)` row.
///
/// The cell, space group, occupancies and displacement parameters come from the
/// structure itself; the model is taken to be the asymmetric unit.
#[pyfunction]
fn structure_factors<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    hkl: &Bound<'py, PyArray2<i32>>,
) -> PyResult<Bound<'py, PyArray1<Complex64>>> {
    let hkl = hkl.readonly();
    let shape = hkl.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("hkl must have shape (n, 3)"));
    }
    let rows: Vec<[i32; 3]> = hkl
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    let structure = structure.inner.clone();
    let values = py
        .detach(move || molframe::crystal::structure_factors(structure.engine(), &rows))
        .map_err(|diagnostic| PyValueError::new_err(diagnostic.to_string()))?;
    let values: Vec<Complex64> = values
        .into_iter()
        .map(|value| Complex64::new(value.re, value.im))
        .collect();
    Ok(values.to_pyarray(py))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(structure_factors, module)?)?;
    module.add_class::<PyUnitCell>()?;
    module.add_class::<PySpaceGroup>()?;
    module.add_class::<PyReflectionSymmetry>()
}
