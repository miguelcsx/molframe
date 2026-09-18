//! Typed Python values for Tripos MOL2 records.

use super::molecules::PyMolecule;
use ::molframe::chem as molframe;
use pyo3::prelude::*;

#[pyclass(name = "Mol2AtomMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMol2AtomMetadata(pub(crate) molframe::Mol2AtomMetadata);

#[pymethods]
impl PyMol2AtomMetadata {
    #[new]
    #[pyo3(signature = (id, name, atom_type, *, substructure_id=None, substructure_name=None, charge=None, status_bits=None))]
    fn new(
        id: usize,
        name: String,
        atom_type: String,
        substructure_id: Option<usize>,
        substructure_name: Option<String>,
        charge: Option<f64>,
        status_bits: Option<String>,
    ) -> PyResult<Self> {
        if id == 0
            || name.is_empty()
            || atom_type.is_empty()
            || name.contains(char::is_whitespace)
            || atom_type.contains(char::is_whitespace)
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MOL2 atom metadata has an invalid identifier or token",
            ));
        }
        if charge.is_some_and(|value| !value.is_finite()) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MOL2 atom charge must be finite",
            ));
        }
        Ok(Self(molframe::Mol2AtomMetadata {
            id,
            name: name.into_boxed_str(),
            atom_type: atom_type.into_boxed_str(),
            substructure_id,
            substructure_name: substructure_name.map(Into::into),
            charge,
            status_bits: status_bits.map(Into::into),
        }))
    }

    #[getter]
    const fn id(&self) -> usize {
        self.0.id
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn atom_type(&self) -> &str {
        &self.0.atom_type
    }

    #[getter]
    const fn substructure_id(&self) -> Option<usize> {
        self.0.substructure_id
    }

    #[getter]
    fn substructure_name(&self) -> Option<&str> {
        self.0.substructure_name.as_deref()
    }

    #[getter]
    const fn charge(&self) -> Option<f64> {
        self.0.charge
    }

    #[getter]
    fn status_bits(&self) -> Option<&str> {
        self.0.status_bits.as_deref()
    }
}

#[pyclass(name = "Mol2BondMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMol2BondMetadata(pub(crate) molframe::Mol2BondMetadata);

#[pymethods]
impl PyMol2BondMetadata {
    #[new]
    #[pyo3(signature = (id, bond_type, *, status_bits=None))]
    fn new(id: usize, bond_type: String, status_bits: Option<String>) -> PyResult<Self> {
        if id == 0 || bond_type.is_empty() || bond_type.contains(char::is_whitespace) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MOL2 bond metadata has an invalid identifier or token",
            ));
        }
        Ok(Self(molframe::Mol2BondMetadata {
            id,
            bond_type: bond_type.into_boxed_str(),
            status_bits: status_bits.map(Into::into),
        }))
    }

    #[getter]
    const fn id(&self) -> usize {
        self.0.id
    }

    #[getter]
    fn bond_type(&self) -> &str {
        &self.0.bond_type
    }

    #[getter]
    fn status_bits(&self) -> Option<&str> {
        self.0.status_bits.as_deref()
    }
}

#[pyclass(name = "Mol2Section", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMol2Section(pub(crate) molframe::Mol2Section);

#[pymethods]
impl PyMol2Section {
    #[new]
    fn new(name: String, lines: Vec<String>) -> PyResult<Self> {
        if name.is_empty() || name.contains(char::is_whitespace) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MOL2 section name must be a non-empty token",
            ));
        }
        Ok(Self(molframe::Mol2Section {
            name: name.into_boxed_str(),
            lines: lines.into_iter().map(Into::into).collect(),
        }))
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn lines(&self) -> Vec<String> {
        self.0.lines.iter().map(ToString::to_string).collect()
    }
}

#[pyclass(name = "Mol2Record", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMol2Record(pub(crate) molframe::Mol2Record);

#[pymethods]
impl PyMol2Record {
    #[new]
    #[pyo3(signature = (name, molecule_type, charge_type, additional_counts, status_bits, comment, molecule, atom_metadata, bond_metadata, extra_sections))]
    fn new(
        name: String,
        molecule_type: String,
        charge_type: String,
        additional_counts: Vec<usize>,
        status_bits: Option<String>,
        comment: Option<String>,
        molecule: &PyMolecule,
        atom_metadata: Vec<PyMol2AtomMetadata>,
        bond_metadata: Vec<PyMol2BondMetadata>,
        extra_sections: Vec<PyMol2Section>,
    ) -> PyResult<Self> {
        if atom_metadata.len() != molecule.0.atoms.len()
            || bond_metadata.len() != molecule.0.bonds.len()
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MOL2 metadata length does not match the molecular graph",
            ));
        }
        Ok(Self(molframe::Mol2Record {
            name: name.into_boxed_str(),
            molecule_type: molecule_type.into_boxed_str(),
            charge_type: charge_type.into_boxed_str(),
            additional_counts,
            status_bits: status_bits.map(Into::into),
            comment: comment.map(Into::into),
            molecule: molecule.0.clone(),
            atom_metadata: atom_metadata.into_iter().map(|item| item.0).collect(),
            bond_metadata: bond_metadata.into_iter().map(|item| item.0).collect(),
            extra_sections: extra_sections.into_iter().map(|item| item.0).collect(),
        }))
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn molecule_type(&self) -> &str {
        &self.0.molecule_type
    }

    #[getter]
    fn charge_type(&self) -> &str {
        &self.0.charge_type
    }

    #[getter]
    fn additional_counts(&self) -> Vec<usize> {
        self.0.additional_counts.clone()
    }

    #[getter]
    fn status_bits(&self) -> Option<&str> {
        self.0.status_bits.as_deref()
    }

    #[getter]
    fn comment(&self) -> Option<&str> {
        self.0.comment.as_deref()
    }

    #[getter]
    fn molecule(&self) -> PyMolecule {
        PyMolecule(self.0.molecule.clone())
    }

    #[getter]
    fn atom_metadata(&self) -> Vec<PyMol2AtomMetadata> {
        self.0
            .atom_metadata
            .iter()
            .cloned()
            .map(PyMol2AtomMetadata)
            .collect()
    }

    #[getter]
    fn bond_metadata(&self) -> Vec<PyMol2BondMetadata> {
        self.0
            .bond_metadata
            .iter()
            .cloned()
            .map(PyMol2BondMetadata)
            .collect()
    }

    #[getter]
    fn extra_sections(&self) -> Vec<PyMol2Section> {
        self.0
            .extra_sections
            .iter()
            .cloned()
            .map(PyMol2Section)
            .collect()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMol2AtomMetadata>()?;
    module.add_class::<PyMol2BondMetadata>()?;
    module.add_class::<PyMol2Section>()?;
    module.add_class::<PyMol2Record>()?;
    Ok(())
}
