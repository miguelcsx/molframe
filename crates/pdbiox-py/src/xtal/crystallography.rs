//! Vectorized unit-cell transforms and the native Hall-setting catalogue.

use crate::structure::PyStructure;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "UnitCell", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyUnitCell {
    pub(crate) cell: pdbiox::UnitCell,
    transform: pdbiox::CellTransform,
}

#[pymethods]
impl PyUnitCell {
    #[new]
    fn new(lengths: [f64; 3], angles: [f64; 3]) -> PyResult<Self> {
        Self::build(pdbiox::UnitCell { lengths, angles })
    }

    #[getter]
    const fn lengths(&self) -> [f64; 3] {
        self.cell.lengths
    }
    #[getter]
    const fn angles(&self) -> [f64; 3] {
        self.cell.angles
    }

    fn is_placeholder(&self) -> bool {
        self.cell.is_placeholder()
    }
    #[getter]
    const fn fractional_to_cartesian_matrix(&self) -> [[f64; 3]; 3] {
        *self.transform.forward_matrix()
    }
    #[getter]
    const fn cartesian_to_fractional_matrix(&self) -> [[f64; 3]; 3] {
        *self.transform.inverse_matrix()
    }

    fn to_cartesian<'py>(
        &self,
        py: Python<'py>,
        coordinates: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let coordinates = points(coordinates)?;
        let transform = self.transform;
        let output = py.detach(move || {
            coordinates
                .into_iter()
                .map(|point| transform.to_cartesian(point))
                .collect::<Vec<_>>()
        });
        array(py, output)
    }

    fn to_fractional<'py>(
        &self,
        py: Python<'py>,
        coordinates: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let coordinates = points(coordinates)?;
        let transform = self.transform;
        let output = py.detach(move || {
            coordinates
                .into_iter()
                .map(|point| transform.to_fractional(point))
                .collect::<Vec<_>>()
        });
        array(py, output)
    }
}

impl PyUnitCell {
    pub(crate) fn from_native(cell: pdbiox::UnitCell) -> PyResult<Self> {
        Self::build(cell)
    }

    fn build(cell: pdbiox::UnitCell) -> PyResult<Self> {
        pdbiox::CellTransform::new(&cell)
            .map(|transform| Self { cell, transform })
            .map_err(value_error)
    }
}

#[pyclass(name = "CellTransform", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCellTransform(pdbiox::CellTransform);

#[pymethods]
impl PyCellTransform {
    #[new]
    fn new(cell: &PyUnitCell) -> PyResult<Self> {
        pdbiox::CellTransform::new(&cell.cell)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn fractional_to_cartesian_matrix(&self) -> [[f64; 3]; 3] {
        *self.0.forward_matrix()
    }

    #[getter]
    fn cartesian_to_fractional_matrix(&self) -> [[f64; 3]; 3] {
        *self.0.inverse_matrix()
    }

    #[pyo3(name = "to_cartesian")]
    fn cartesian(&self, fractional: [f64; 3]) -> [f64; 3] {
        self.0.to_cartesian(fractional)
    }

    #[pyo3(name = "to_fractional")]
    fn fractional(&self, cartesian: [f64; 3]) -> [f64; 3] {
        self.0.to_fractional(cartesian)
    }
}

#[pyclass(name = "SymmetryOperation", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySymmetryOperation {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    rotation: [[i32; 3]; 3],
    #[pyo3(get)]
    translation: [(i32, u32); 3],
}

#[pyclass(name = "SpaceGroup", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySpaceGroup {
    #[pyo3(get)]
    hall_number: u16,
    #[pyo3(get)]
    international_number: u16,
    #[pyo3(get)]
    international_short: String,
    #[pyo3(get)]
    international_full: String,
    #[pyo3(get)]
    hall_symbol: String,
    #[pyo3(get)]
    choice: String,
    #[pyo3(get)]
    operations: Vec<PySymmetryOperation>,
}

#[pyfunction]
pub(crate) fn space_group_by_number(hall_number: u16) -> PyResult<PySpaceGroup> {
    pdbiox::space_group_setting(hall_number)
        .map(PySpaceGroup::from)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn space_group_by_symbol(hall_symbol: &str) -> PyResult<PySpaceGroup> {
    pdbiox::space_group_by_hall(hall_symbol)
        .map(PySpaceGroup::from)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn space_group_setting(hall_number: u16) -> PyResult<PySpaceGroup> {
    space_group_by_number(hall_number)
}

#[pyfunction]
pub(crate) fn space_group_by_hall(hall_symbol: &str) -> PyResult<PySpaceGroup> {
    space_group_by_symbol(hall_symbol)
}

#[pyfunction]
pub(crate) fn space_group_settings(international_number: u16) -> PyResult<Vec<PySpaceGroup>> {
    pdbiox::space_group_settings(international_number)
        .map(|values| values.into_iter().map(PySpaceGroup::from).collect())
        .map_err(value_error)
}

#[pymethods]
impl PyStructure {
    #[getter]
    fn cell(&self) -> PyResult<Option<PyUnitCell>> {
        self.structure()
            .data()
            .cell
            .map(PyUnitCell::build)
            .transpose()
    }
}

impl From<&pdbiox::SpaceGroupSetting> for PySpaceGroup {
    fn from(value: &pdbiox::SpaceGroupSetting) -> Self {
        Self {
            hall_number: value.hall_number,
            international_number: value.international_number,
            international_short: value.international_short.to_string(),
            international_full: value.international_full.to_string(),
            hall_symbol: value.hall_symbol.to_string(),
            choice: value.choice.to_string(),
            operations: value
                .operations
                .iter()
                .map(PySymmetryOperation::from)
                .collect(),
        }
    }
}

impl From<&pdbiox::SymmetryOperation> for PySymmetryOperation {
    fn from(value: &pdbiox::SymmetryOperation) -> Self {
        Self {
            id: value.id.to_string(),
            rotation: value.rotation,
            translation: value
                .translation
                .map(|part| (part.numerator(), part.denominator())),
        }
    }
}

fn points(values: PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    let shape = values.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("coordinates must have shape (n, 3)"));
    }
    let output = values
        .as_slice()?
        .chunks_exact(3)
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    drop(values);
    Ok(output)
}

fn array(py: Python<'_>, values: Vec<[f64; 3]>) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let rows = values.len();
    let flat = values.into_iter().flatten().collect();
    Array2::from_shape_vec((rows, 3), flat)
        .map(|values| values.into_pyarray(py))
        .map_err(value_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::xtal_maps::register(module)?;
    crate::xtal_restraints::register(module)?;
    module.add_class::<PyUnitCell>()?;
    module.add_class::<PyCellTransform>()?;
    module.add_class::<PySymmetryOperation>()?;
    module.add_class::<PySpaceGroup>()?;
    module.add("SpaceGroupSetting", module.getattr("SpaceGroup")?)?;
    module.add("ASSEMBLIES_EXTENSION", pdbiox::xtal::ASSEMBLIES_EXTENSION)?;
    module.add(
        "INSTANCE_ID_ANNOTATION",
        pdbiox::xtal::INSTANCE_ID_ANNOTATION,
    )?;
    module.add("NCS_EXTENSION", pdbiox::xtal::NCS_EXTENSION)?;
    module.add("SYMMETRY_EXTENSION", pdbiox::xtal::SYMMETRY_EXTENSION)?;
    module.add_function(wrap_pyfunction!(space_group_by_number, module)?)?;
    module.add_function(wrap_pyfunction!(space_group_by_symbol, module)?)?;
    module.add_function(wrap_pyfunction!(space_group_by_hall, module)?)?;
    module.add_function(wrap_pyfunction!(space_group_setting, module)?)?;
    module.add_function(wrap_pyfunction!(space_group_settings, module)?)?;
    crate::xtal::register(module)?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
