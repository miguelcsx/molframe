//! Python projections for monomer restraint libraries and structure-factor CIF.

use crate::cif_document::PyCifDocument;
use crate::xtal::PyReflectionTable;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "BondRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondRestraint {
    #[pyo3(get)]
    atoms: [String; 2],
    #[pyo3(get)]
    target: f64,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    kind: Option<String>,
}

#[pyclass(name = "AngleRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAngleRestraint {
    #[pyo3(get)]
    atoms: [String; 3],
    #[pyo3(get)]
    target: f64,
    #[pyo3(get)]
    sigma: f64,
}

#[pyclass(name = "TorsionRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTorsionRestraint {
    #[pyo3(get)]
    id: Option<String>,
    #[pyo3(get)]
    atoms: [String; 4],
    #[pyo3(get)]
    target: f64,
    #[pyo3(get)]
    sigma: f64,
    #[pyo3(get)]
    period: u32,
}

#[pyclass(name = "PlaneAtomRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlaneAtomRestraint {
    #[pyo3(get)]
    atom: String,
    #[pyo3(get)]
    sigma: f64,
}

#[pyclass(name = "PlaneRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlaneRestraint {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    atoms: Vec<PyPlaneAtomRestraint>,
}

#[pyclass(name = "ChiralVolumeSign", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyChiralVolumeSign {
    Positive,
    Negative,
    Both,
}

#[pyclass(name = "ChiralRestraint", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChiralRestraint {
    #[pyo3(get)]
    id: Option<String>,
    #[pyo3(get)]
    atoms: [String; 4],
    #[pyo3(get)]
    sign: PyChiralVolumeSign,
}

#[pyclass(name = "MonomerRestraints", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMonomerRestraints {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    atoms: Vec<String>,
    #[pyo3(get)]
    bonds: Vec<PyBondRestraint>,
    #[pyo3(get)]
    angles: Vec<PyAngleRestraint>,
    #[pyo3(get)]
    torsions: Vec<PyTorsionRestraint>,
    #[pyo3(get)]
    planes: Vec<PyPlaneRestraint>,
    #[pyo3(get)]
    chirals: Vec<PyChiralRestraint>,
}

#[pyclass(name = "MonomerLibrary", frozen, skip_from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyMonomerLibrary(pub(crate) pdbiox::xtal::MonomerLibrary);

#[pymethods]
impl PyMonomerLibrary {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn get(&self, id: &str) -> Option<PyMonomerRestraints> {
        self.0.get(id).map(Into::into)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn items(&self) -> Vec<(String, PyMonomerRestraints)> {
        self.0
            .iter()
            .map(|(id, restraints)| (id.to_owned(), restraints.into()))
            .collect()
    }
}

#[pyfunction]
pub(crate) fn lower_monomer_library(
    py: Python<'_>,
    document: &PyCifDocument,
) -> PyResult<PyMonomerLibrary> {
    py.detach(move || -> PyResult<PyMonomerLibrary> {
        pdbiox::xtal::lower_monomer_library(&document.inner)
            .map(PyMonomerLibrary)
            .map_err(value_error)
    })
}

#[pyfunction]
pub(crate) fn read_monomer_library(
    py: Python<'_>,
    data: &[u8],
) -> PyResult<(PyMonomerLibrary, Vec<String>)> {
    let input = pdbiox::InputBuffer::from_bytes(data.to_vec());
    py.detach(move || pdbiox::xtal::read_monomer_library(&input))
        .map(|(library, findings)| {
            (
                PyMonomerLibrary(library),
                findings
                    .into_iter()
                    .map(|finding| finding.to_string())
                    .collect(),
            )
        })
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn lower_structure_factor_cif(
    py: Python<'_>,
    document: &PyCifDocument,
) -> PyResult<PyReflectionTable> {
    py.detach(|| {
        pdbiox::xtal::lower_structure_factor_cif(&document.inner)
            .map(PyReflectionTable)
            .map_err(value_error)
    })
}

#[pyfunction]
pub(crate) fn write_structure_factor_cif(
    py: Python<'_>,
    table: &PyReflectionTable,
) -> PyResult<String> {
    py.detach(|| pdbiox::xtal::write_structure_factor_cif(&table.0).map_err(value_error))
}

impl From<&pdbiox::xtal::MonomerRestraints> for PyMonomerRestraints {
    fn from(value: &pdbiox::xtal::MonomerRestraints) -> Self {
        Self {
            id: value.id.to_string(),
            atoms: value.atoms.iter().map(ToString::to_string).collect(),
            bonds: value.bonds.iter().map(Into::into).collect(),
            angles: value.angles.iter().map(Into::into).collect(),
            torsions: value.torsions.iter().map(Into::into).collect(),
            planes: value.planes.iter().map(Into::into).collect(),
            chirals: value.chirals.iter().map(Into::into).collect(),
        }
    }
}

impl From<&pdbiox::xtal::BondRestraint> for PyBondRestraint {
    fn from(value: &pdbiox::xtal::BondRestraint) -> Self {
        Self {
            atoms: value.atoms.clone().map(|atom| atom.to_string()),
            target: value.target,
            sigma: value.sigma,
            kind: value.kind.as_deref().map(str::to_owned),
        }
    }
}

impl From<&pdbiox::xtal::AngleRestraint> for PyAngleRestraint {
    fn from(value: &pdbiox::xtal::AngleRestraint) -> Self {
        Self {
            atoms: value.atoms.clone().map(|atom| atom.to_string()),
            target: value.target,
            sigma: value.sigma,
        }
    }
}

impl From<&pdbiox::xtal::TorsionRestraint> for PyTorsionRestraint {
    fn from(value: &pdbiox::xtal::TorsionRestraint) -> Self {
        Self {
            id: value.id.as_deref().map(str::to_owned),
            atoms: value.atoms.clone().map(|atom| atom.to_string()),
            target: value.target,
            sigma: value.sigma,
            period: value.period,
        }
    }
}

impl From<&pdbiox::xtal::PlaneAtomRestraint> for PyPlaneAtomRestraint {
    fn from(value: &pdbiox::xtal::PlaneAtomRestraint) -> Self {
        Self {
            atom: value.atom.to_string(),
            sigma: value.sigma,
        }
    }
}

impl From<&pdbiox::xtal::PlaneRestraint> for PyPlaneRestraint {
    fn from(value: &pdbiox::xtal::PlaneRestraint) -> Self {
        Self {
            id: value.id.to_string(),
            atoms: value.atoms.iter().map(Into::into).collect(),
        }
    }
}

impl From<pdbiox::xtal::ChiralVolumeSign> for PyChiralVolumeSign {
    fn from(value: pdbiox::xtal::ChiralVolumeSign) -> Self {
        match value {
            pdbiox::xtal::ChiralVolumeSign::Positive => Self::Positive,
            pdbiox::xtal::ChiralVolumeSign::Negative => Self::Negative,
            pdbiox::xtal::ChiralVolumeSign::Both => Self::Both,
        }
    }
}

impl From<&pdbiox::xtal::ChiralRestraint> for PyChiralRestraint {
    fn from(value: &pdbiox::xtal::ChiralRestraint) -> Self {
        Self {
            id: value.id.as_deref().map(str::to_owned),
            atoms: value.atoms.clone().map(|atom| atom.to_string()),
            sign: value.sign.into(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBondRestraint>()?;
    module.add_class::<PyAngleRestraint>()?;
    module.add_class::<PyTorsionRestraint>()?;
    module.add_class::<PyPlaneAtomRestraint>()?;
    module.add_class::<PyPlaneRestraint>()?;
    module.add_class::<PyChiralVolumeSign>()?;
    module.add_class::<PyChiralRestraint>()?;
    module.add_class::<PyMonomerRestraints>()?;
    module.add_class::<PyMonomerLibrary>()?;
    module.add_function(wrap_pyfunction!(lower_monomer_library, module)?)?;
    module.add_function(wrap_pyfunction!(read_monomer_library, module)?)?;
    module.add_function(wrap_pyfunction!(lower_structure_factor_cif, module)?)?;
    module.add_function(wrap_pyfunction!(write_structure_factor_cif, module)?)?;
    Ok(())
}
