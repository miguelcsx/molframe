//! Mechanical Python handle over the Rust bond table.

use crate::errors::index_error;
use crate::index::{PyAtomIndex, PyBondIndex, normalise_index};
use crate::structure::PyStructure;
use pdbiox::Structure;
use pyo3::prelude::*;

#[pyclass(name = "BondOrder", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyBondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Polymeric,
    Unknown,
}

#[pyclass(name = "BondProvenance", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyBondProvenance {
    File,
    ChemicalComponentDictionary,
    InferredDistance,
    User,
}

#[pyclass(name = "BondRecord", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBondRecord {
    #[pyo3(get)]
    atom_a: u32,
    #[pyo3(get)]
    atom_b: u32,
    #[pyo3(get)]
    order: PyBondOrder,
    #[pyo3(get)]
    provenance: PyBondProvenance,
}

impl PyBondRecord {
    fn native(self) -> pdbiox::BondRecord {
        pdbiox::BondRecord {
            atom_a: pdbiox::AtomIndex::new(self.atom_a),
            atom_b: pdbiox::AtomIndex::new(self.atom_b),
            order: self.order.into(),
            provenance: self.provenance.into(),
        }
    }
}

#[pyclass(name = "BondTable", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondTable(pub(crate) pdbiox::BondTable);

#[pyclass(name = "BondTableBuilder", unsendable, skip_from_py_object)]
pub(crate) struct PyBondTableBuilder(pdbiox::BondTableBuilder);

#[pyclass(name = "BondAdjacency", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondAdjacency(pub(crate) pdbiox::BondAdjacency);

#[pyclass(name = "Bonds", frozen, skip_from_py_object)]
pub(crate) struct PyBonds {
    inner: Structure,
}

impl PyBonds {
    pub(crate) const fn new(inner: Structure) -> Self {
        Self { inner }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }
}

#[pymethods]
impl PyBondRecord {
    #[new]
    #[pyo3(signature = (atom_a, atom_b, order=PyBondOrder::Unknown, provenance=PyBondProvenance::User))]
    fn new(
        atom_a: PyAtomIndex,
        atom_b: PyAtomIndex,
        order: PyBondOrder,
        provenance: PyBondProvenance,
    ) -> Self {
        Self {
            atom_a: atom_a.0.get(),
            atom_b: atom_b.0.get(),
            order,
            provenance,
        }
    }

    #[getter]
    fn atom_a_index(&self) -> PyAtomIndex {
        PyAtomIndex(pdbiox::AtomIndex::new(self.atom_a))
    }

    #[getter]
    fn atom_b_index(&self) -> PyAtomIndex {
        PyAtomIndex(pdbiox::AtomIndex::new(self.atom_b))
    }
}

#[pymethods]
impl PyBondTable {
    #[new]
    fn new() -> Self {
        Self(pdbiox::BondTable::default())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[getter]
    fn available(&self) -> bool {
        self.0.is_available()
    }

    fn is_available(&self) -> bool {
        self.0.is_available()
    }

    fn get(&self, bond: PyBondIndex) -> Option<PyBondRecord> {
        self.0.get(bond.0).map(PyBondRecord::from)
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<PyBondRecord> {
        let Some(position) = normalise_index(index, self.0.len()) else {
            return Err(index_error(py, index));
        };
        let position = u32::try_from(position).map_err(|_| index_error(py, index))?;
        self.0
            .get(pdbiox::BondIndex::new(position))
            .map(PyBondRecord::from)
            .ok_or_else(|| index_error(py, index))
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> Vec<PyBondRecord> {
        self.0.iter().map(PyBondRecord::from).collect()
    }

    fn adjacency(&self, atom_count: u32) -> PyBondAdjacency {
        PyBondAdjacency(self.0.adjacency(atom_count).clone())
    }
}

#[pymethods]
impl PyBondTableBuilder {
    #[new]
    fn new() -> Self {
        Self(pdbiox::BondTableBuilder::new())
    }

    fn push(&mut self, record: PyBondRecord) {
        self.0.push(record.native());
    }

    fn finish(&mut self) -> PyBondTable {
        PyBondTable(std::mem::take(&mut self.0).finish())
    }

    fn finish_with_availability(&mut self, available: bool) -> PyBondTable {
        PyBondTable(std::mem::take(&mut self.0).finish_with_availability(available))
    }
}

#[pymethods]
impl PyBondAdjacency {
    fn neighbours(&self, atom: PyAtomIndex) -> Vec<PyAtomIndex> {
        self.0
            .neighbours(atom.0)
            .iter()
            .copied()
            .map(PyAtomIndex)
            .collect()
    }
}

#[pymethods]
impl PyStructure {
    #[getter]
    fn bonds(&self) -> PyBonds {
        PyBonds::new(self.structure().clone())
    }
}

#[pymethods]
impl PyBonds {
    fn __len__(&self) -> usize {
        self.inner.data().bonds.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<PyBondRecord> {
        let Some(position) = normalise_index(index, self.__len__()) else {
            return Err(index_error(py, index));
        };
        let position = u32::try_from(position).map_err(|_| index_error(py, index))?;
        self.inner
            .data()
            .bonds
            .get(pdbiox::BondIndex::new(position))
            .map(PyBondRecord::from)
            .ok_or_else(|| index_error(py, index))
    }

    #[getter]
    fn available(&self) -> bool {
        self.inner.data().bonds.is_available()
    }
}

impl From<pdbiox::BondRecord> for PyBondRecord {
    fn from(value: pdbiox::BondRecord) -> Self {
        Self {
            atom_a: value.atom_a.get(),
            atom_b: value.atom_b.get(),
            order: value.order.into(),
            provenance: value.provenance.into(),
        }
    }
}

impl From<pdbiox::BondOrder> for PyBondOrder {
    fn from(value: pdbiox::BondOrder) -> Self {
        match value {
            pdbiox::BondOrder::Single => Self::Single,
            pdbiox::BondOrder::Double => Self::Double,
            pdbiox::BondOrder::Triple => Self::Triple,
            pdbiox::BondOrder::Quadruple => Self::Quadruple,
            pdbiox::BondOrder::Aromatic => Self::Aromatic,
            pdbiox::BondOrder::Polymeric => Self::Polymeric,
            pdbiox::BondOrder::Unknown => Self::Unknown,
        }
    }
}

impl From<PyBondOrder> for pdbiox::BondOrder {
    fn from(value: PyBondOrder) -> Self {
        match value {
            PyBondOrder::Single => Self::Single,
            PyBondOrder::Double => Self::Double,
            PyBondOrder::Triple => Self::Triple,
            PyBondOrder::Quadruple => Self::Quadruple,
            PyBondOrder::Aromatic => Self::Aromatic,
            PyBondOrder::Polymeric => Self::Polymeric,
            PyBondOrder::Unknown => Self::Unknown,
        }
    }
}

impl From<pdbiox::BondProvenance> for PyBondProvenance {
    fn from(value: pdbiox::BondProvenance) -> Self {
        match value {
            pdbiox::BondProvenance::File => Self::File,
            pdbiox::BondProvenance::ChemicalComponentDictionary => {
                Self::ChemicalComponentDictionary
            }
            pdbiox::BondProvenance::InferredDistance => Self::InferredDistance,
            pdbiox::BondProvenance::User => Self::User,
        }
    }
}

impl From<PyBondProvenance> for pdbiox::BondProvenance {
    fn from(value: PyBondProvenance) -> Self {
        match value {
            PyBondProvenance::File => Self::File,
            PyBondProvenance::ChemicalComponentDictionary => Self::ChemicalComponentDictionary,
            PyBondProvenance::InferredDistance => Self::InferredDistance,
            PyBondProvenance::User => Self::User,
        }
    }
}
