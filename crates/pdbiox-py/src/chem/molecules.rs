//! Typed Python values for MOL, SDF and MOL2 records.

use crate::chemistry::PyElement;
use ::pdbiox::chem as pdbiox;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, MolError, PyValueError);
create_exception!(_native, Mol2Error, PyValueError);

#[pyclass(name = "MolAtom", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMolAtom(pub(crate) pdbiox::MolAtom);

#[pymethods]
impl PyMolAtom {
    #[new]
    fn new(element: &PyElement, position: [f32; 3]) -> PyResult<Self> {
        if !position.iter().all(|value| value.is_finite()) {
            return Err(PyValueError::new_err("atom position must be finite"));
        }
        Ok(Self(pdbiox::MolAtom {
            element: element.0,
            position,
        }))
    }

    #[getter]
    fn element(&self) -> PyElement {
        PyElement(self.0.element)
    }

    #[getter]
    const fn position(&self) -> [f32; 3] {
        self.0.position
    }
}

#[pyclass(name = "MolBond", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMolBond(pub(crate) pdbiox::MolBond);

#[pymethods]
impl PyMolBond {
    #[new]
    fn new(first: usize, second: usize, order: u8) -> PyResult<Self> {
        if !(1..=4).contains(&order) {
            return Err(PyValueError::new_err("MOL bond order must be in 1..=4"));
        }
        Ok(Self(pdbiox::MolBond {
            first,
            second,
            order,
        }))
    }

    #[getter]
    const fn first(&self) -> usize {
        self.0.first
    }

    #[getter]
    const fn second(&self) -> usize {
        self.0.second
    }

    #[getter]
    const fn order(&self) -> u8 {
        self.0.order
    }
}

#[pyclass(name = "Molecule", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMolecule(pub(crate) pdbiox::Molecule);

#[pymethods]
impl PyMolecule {
    #[new]
    fn new(atoms: Vec<PyMolAtom>, bonds: Vec<PyMolBond>) -> PyResult<Self> {
        let atoms = atoms.into_iter().map(|atom| atom.0).collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(|bond| bond.0).collect::<Vec<_>>();
        validate_bonds(atoms.len(), &bonds)?;
        Ok(Self(pdbiox::Molecule { atoms, bonds }))
    }

    #[getter]
    fn atoms(&self) -> Vec<PyMolAtom> {
        self.0.atoms.iter().copied().map(PyMolAtom).collect()
    }

    #[getter]
    fn bonds(&self) -> Vec<PyMolBond> {
        self.0.bonds.iter().copied().map(PyMolBond).collect()
    }
}

#[pyclass(name = "MolVersion", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMolVersion {
    V2000,
    V3000,
}

impl From<PyMolVersion> for pdbiox::MolVersion {
    fn from(value: PyMolVersion) -> Self {
        match value {
            PyMolVersion::V2000 => Self::V2000,
            PyMolVersion::V3000 => Self::V3000,
        }
    }
}

impl From<pdbiox::MolVersion> for PyMolVersion {
    fn from(value: pdbiox::MolVersion) -> Self {
        match value {
            pdbiox::MolVersion::V2000 => Self::V2000,
            pdbiox::MolVersion::V3000 => Self::V3000,
        }
    }
}

#[pyclass(name = "MolAtomMetadata", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMolAtomMetadata(pub(crate) pdbiox::MolAtomMetadata);

#[pymethods]
impl PyMolAtomMetadata {
    #[new]
    #[pyo3(signature = (*, formal_charge=None, isotope=None, stereo_parity=None))]
    fn new(formal_charge: Option<i8>, isotope: Option<u16>, stereo_parity: Option<u8>) -> Self {
        Self(pdbiox::MolAtomMetadata {
            formal_charge,
            isotope,
            stereo_parity,
        })
    }

    #[getter]
    const fn formal_charge(&self) -> Option<i8> {
        self.0.formal_charge
    }

    #[getter]
    const fn isotope(&self) -> Option<u16> {
        self.0.isotope
    }

    #[getter]
    const fn stereo_parity(&self) -> Option<u8> {
        self.0.stereo_parity
    }
}

#[pyclass(name = "MolBondMetadata", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMolBondMetadata(pub(crate) pdbiox::MolBondMetadata);

#[pymethods]
impl PyMolBondMetadata {
    #[new]
    #[pyo3(signature = (*, stereo=None))]
    fn new(stereo: Option<u8>) -> Self {
        Self(pdbiox::MolBondMetadata { stereo })
    }

    #[getter]
    const fn stereo(&self) -> Option<u8> {
        self.0.stereo
    }
}

#[pyclass(name = "SdfProperty", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySdfProperty(pub(crate) pdbiox::SdfProperty);

#[pymethods]
impl PySdfProperty {
    #[new]
    fn new(name: String, value: String) -> PyResult<Self> {
        if name.is_empty() || name.contains(['<', '>', '\n', '\r']) {
            return Err(PyValueError::new_err("invalid SDF property name"));
        }
        Ok(Self(pdbiox::SdfProperty {
            name: name.into_boxed_str(),
            value: value.into_boxed_str(),
        }))
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn value(&self) -> &str {
        &self.0.value
    }
}

#[pyclass(name = "MolRecord", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMolRecord(pub(crate) pdbiox::MolRecord);

#[pymethods]
impl PyMolRecord {
    #[new]
    #[pyo3(signature = (name, program, comment, version, molecule, atom_metadata, bond_metadata, properties))]
    fn new(
        name: String,
        program: String,
        comment: String,
        version: PyMolVersion,
        molecule: &PyMolecule,
        atom_metadata: Vec<PyMolAtomMetadata>,
        bond_metadata: Vec<PyMolBondMetadata>,
        properties: Vec<PySdfProperty>,
    ) -> PyResult<Self> {
        if name.contains(['\n', '\r'])
            || program.contains(['\n', '\r'])
            || comment.contains(['\n', '\r'])
            || atom_metadata.len() != molecule.0.atoms.len()
            || bond_metadata.len() != molecule.0.bonds.len()
        {
            return Err(PyValueError::new_err("invalid MOL record fields"));
        }
        Ok(Self(pdbiox::MolRecord {
            name: name.into_boxed_str(),
            program: program.into_boxed_str(),
            comment: comment.into_boxed_str(),
            version: version.into(),
            molecule: molecule.0.clone(),
            atom_metadata: atom_metadata.into_iter().map(|item| item.0).collect(),
            bond_metadata: bond_metadata.into_iter().map(|item| item.0).collect(),
            properties: properties.into_iter().map(|item| item.0).collect(),
        }))
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn program(&self) -> &str {
        &self.0.program
    }

    #[getter]
    fn comment(&self) -> &str {
        &self.0.comment
    }

    #[getter]
    fn version(&self) -> PyMolVersion {
        self.0.version.into()
    }

    #[getter]
    fn molecule(&self) -> PyMolecule {
        PyMolecule(self.0.molecule.clone())
    }

    #[getter]
    fn atom_metadata(&self) -> Vec<PyMolAtomMetadata> {
        self.0
            .atom_metadata
            .iter()
            .copied()
            .map(PyMolAtomMetadata)
            .collect()
    }

    #[getter]
    fn bond_metadata(&self) -> Vec<PyMolBondMetadata> {
        self.0
            .bond_metadata
            .iter()
            .copied()
            .map(PyMolBondMetadata)
            .collect()
    }

    #[getter]
    fn properties(&self) -> Vec<PySdfProperty> {
        self.0
            .properties
            .iter()
            .cloned()
            .map(PySdfProperty)
            .collect()
    }
}

fn validate_bonds(atom_count: usize, bonds: &[pdbiox::MolBond]) -> PyResult<()> {
    if bonds.iter().any(|bond| {
        bond.first >= atom_count || bond.second >= atom_count || !(1..=4).contains(&bond.order)
    }) {
        return Err(PyValueError::new_err(
            "MOL bond endpoint or order is invalid",
        ));
    }
    Ok(())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("MolError", module.py().get_type::<MolError>())?;
    module.add("Mol2Error", module.py().get_type::<Mol2Error>())?;
    module.add_class::<PyMolAtom>()?;
    module.add_class::<PyMolBond>()?;
    module.add_class::<PyMolecule>()?;
    module.add_class::<PyMolVersion>()?;
    module.add_class::<PyMolAtomMetadata>()?;
    module.add_class::<PyMolBondMetadata>()?;
    module.add_class::<PySdfProperty>()?;
    module.add_class::<PyMolRecord>()?;
    super::molecule2::register(module)?;
    super::molecule_io::register(module)
}
