//! Native small-molecule CIF projections.

use crate::cif_document::PyCifDocument;
use crate::crystallography::PyUnitCell;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, SmallCifError, PyValueError);

#[pyclass(name = "SmallCifAtom", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySmallCifAtom {
    #[pyo3(get)]
    label: String,
    #[pyo3(get)]
    type_symbol: String,
    #[pyo3(get)]
    fractional: Option<[f64; 3]>,
    #[pyo3(get)]
    cartesian: Option<[f64; 3]>,
    #[pyo3(get)]
    occupancy: Option<f64>,
    #[pyo3(get)]
    u_iso: Option<f64>,
    #[pyo3(get)]
    charge: Option<String>,
}

#[pyclass(name = "SmallCifBond", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySmallCifBond {
    #[pyo3(get)]
    first: usize,
    #[pyo3(get)]
    second: usize,
    #[pyo3(get)]
    distance: Option<f64>,
    #[pyo3(get)]
    site_symmetry: Option<String>,
    #[pyo3(get)]
    bond_type: Option<String>,
}

#[pyclass(name = "SmallCifStructure", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySmallCifStructure {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    cell: Option<PyUnitCell>,
    #[pyo3(get)]
    space_group_name: Option<String>,
    #[pyo3(get)]
    space_group_number: Option<i32>,
    #[pyo3(get)]
    symmetry_operations: Vec<String>,
    #[pyo3(get)]
    atoms: Vec<PySmallCifAtom>,
    #[pyo3(get)]
    bonds: Vec<PySmallCifBond>,
}

#[pyclass(name = "SmallCifDialect", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySmallCifDialect {
    Ddl2,
    Ddl1,
}

#[pyclass(name = "SmallCifOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySmallCifOptions(pdbiox::cif::SmallCifOptions);

#[pymethods]
impl PySmallCifOptions {
    #[new]
    #[pyo3(signature = (dialect=PySmallCifDialect::Ddl2))]
    fn new(dialect: PySmallCifDialect) -> Self {
        match dialect {
            PySmallCifDialect::Ddl2 => Self(pdbiox::cif::SmallCifOptions::ddl2()),
            PySmallCifDialect::Ddl1 => Self(pdbiox::cif::SmallCifOptions::ddl1()),
        }
    }

    #[staticmethod]
    fn ddl2() -> Self {
        Self::new(PySmallCifDialect::Ddl2)
    }

    #[staticmethod]
    fn ddl1() -> Self {
        Self::new(PySmallCifDialect::Ddl1)
    }
}

#[pyfunction(name = "lower_small_cif")]
pub(crate) fn lower_small_cif(document: &PyCifDocument) -> PyResult<PySmallCifStructure> {
    let value = pdbiox::cif::lower_small_cif(&document.inner)
        .map_err(|error| SmallCifError::new_err(error.to_string()))?;
    PySmallCifStructure::try_from(value).map_err(|error| SmallCifError::new_err(error.to_string()))
}

#[pyfunction(name = "lower_small_cif_with_options")]
pub(crate) fn lower_small_cif_with_options(
    document: &PyCifDocument,
    options: &PySmallCifOptions,
) -> PyResult<PySmallCifStructure> {
    let value = pdbiox::cif::lower_small_cif_with_options(&document.inner, options.0)
        .map_err(|error| SmallCifError::new_err(error.to_string()))?;
    PySmallCifStructure::try_from(value).map_err(|error| SmallCifError::new_err(error.to_string()))
}

impl From<pdbiox::cif::SmallCifAtom> for PySmallCifAtom {
    fn from(value: pdbiox::cif::SmallCifAtom) -> Self {
        Self {
            label: value.label.into(),
            type_symbol: value.type_symbol.into(),
            fractional: value.fractional,
            cartesian: value.cartesian,
            occupancy: value.occupancy,
            u_iso: value.u_iso,
            charge: value.charge.map(Into::into),
        }
    }
}

impl From<pdbiox::cif::SmallCifBond> for PySmallCifBond {
    fn from(value: pdbiox::cif::SmallCifBond) -> Self {
        Self {
            first: value.first,
            second: value.second,
            distance: value.distance,
            site_symmetry: value.site_symmetry.map(Into::into),
            bond_type: value.bond_type.map(Into::into),
        }
    }
}

impl TryFrom<pdbiox::cif::SmallCifStructure> for PySmallCifStructure {
    type Error = PyErr;

    fn try_from(value: pdbiox::cif::SmallCifStructure) -> Result<Self, Self::Error> {
        let cell = value.cell.map(PyUnitCell::from_native).transpose()?;
        Ok(Self {
            name: value.name.into(),
            cell,
            space_group_name: value.space_group_name.map(Into::into),
            space_group_number: value.space_group_number,
            symmetry_operations: value
                .symmetry_operations
                .into_iter()
                .map(Into::into)
                .collect(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value.bonds.into_iter().map(Into::into).collect(),
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("SmallCifError", module.py().get_type::<SmallCifError>())?;
    module.add_class::<PySmallCifAtom>()?;
    module.add_class::<PySmallCifBond>()?;
    module.add_class::<PySmallCifStructure>()?;
    module.add_class::<PySmallCifDialect>()?;
    module.add_class::<PySmallCifOptions>()?;
    module.add_function(wrap_pyfunction!(lower_small_cif, module)?)?;
    module.add_function(wrap_pyfunction!(lower_small_cif_with_options, module)?)?;
    Ok(())
}
