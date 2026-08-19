//! Python projections of the facade's offset-based hierarchy tables.

use crate::core_values::{PyOptionalI32, PyOptionalSymbol, PySymbolId};
use crate::index::{PyChainIndex, PyEntityIndex, PyModelIndex, PyResidueIndex};
use crate::metadata::{PyEntityKind, PyPolymerKind};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::ops::Range;

pub(crate) fn table_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn range_pair(range: Range<u32>) -> [u32; 2] {
    [range.start, range.end]
}

pub(crate) fn checked_range(start: u32, end: u32) -> PyResult<Range<u32>> {
    (start <= end)
        .then_some(start..end)
        .ok_or_else(|| PyValueError::new_err("range end must not precede range start"))
}

#[pyclass(name = "ModelTable", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyModelTable(pub(crate) pdbiox::core::topology::ModelTable);

#[pyclass(name = "ChainRecord", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyChainRecord {
    pub(crate) label_asym_id: PySymbolId,
    pub(crate) auth_asym_id: PyOptionalSymbol,
    pub(crate) entity: PyEntityIndex,
    pub(crate) polymer_kind: PyPolymerKind,
}

#[pyclass(name = "ChainTable", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainTable(pub(crate) pdbiox::core::topology::ChainTable);

#[pyclass(name = "ResidueRecord", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyResidueRecord {
    pub(crate) label_comp_id: PySymbolId,
    pub(crate) auth_comp_id: PyOptionalSymbol,
    pub(crate) label_seq_id: PyOptionalI32,
    pub(crate) auth_seq_id: PyOptionalI32,
    pub(crate) ins_code: PyOptionalSymbol,
    pub(crate) het: bool,
}

#[pyclass(name = "ResidueTable", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyResidueTable(pub(crate) pdbiox::core::topology::ResidueTable);

#[pyclass(name = "EntityTable", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEntityTable(pub(crate) pdbiox::core::topology::EntityTable);

#[pymethods]
impl PyModelTable {
    #[new]
    fn new() -> Self {
        Self(pdbiox::core::topology::ModelTable::default())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn push(&mut self, model_num: i32, first_chain: u32, end_chain: u32) -> PyResult<PyModelIndex> {
        self.0
            .push(model_num, checked_range(first_chain, end_chain)?)
            .map(PyModelIndex)
            .map_err(table_error)
    }

    fn model_num(&self, model: PyModelIndex) -> Option<i32> {
        self.0.model_num(model.0)
    }

    fn chains(&self, model: PyModelIndex) -> Option<[u32; 2]> {
        self.0.chains(model.0).map(range_pair)
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> Vec<PyModelIndex> {
        self.0.iter().map(PyModelIndex).collect()
    }
}

#[pymethods]
impl PyChainRecord {
    #[new]
    #[pyo3(signature = (label_asym_id, entity, polymer_kind=PyPolymerKind::None, auth_asym_id=None))]
    fn new(
        label_asym_id: PySymbolId,
        entity: PyEntityIndex,
        polymer_kind: PyPolymerKind,
        auth_asym_id: Option<PySymbolId>,
    ) -> Self {
        Self {
            label_asym_id,
            auth_asym_id: PyOptionalSymbol::new(auth_asym_id),
            entity,
            polymer_kind,
        }
    }

    #[getter]
    fn label_asym_id(&self) -> PySymbolId {
        self.label_asym_id
    }

    #[getter]
    fn auth_asym_id(&self) -> PyOptionalSymbol {
        self.auth_asym_id
    }

    #[getter]
    fn entity(&self) -> PyEntityIndex {
        self.entity
    }

    #[getter]
    fn polymer_kind(&self) -> PyPolymerKind {
        self.polymer_kind
    }
}

#[pymethods]
impl PyChainTable {
    #[new]
    fn new() -> Self {
        Self(pdbiox::core::topology::ChainTable::default())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn push(
        &mut self,
        record: PyChainRecord,
        first_residue: u32,
        end_residue: u32,
    ) -> PyResult<PyChainIndex> {
        let native = pdbiox::core::topology::ChainRecord {
            label_asym_id: record.label_asym_id.0,
            auth_asym_id: record.auth_asym_id.0,
            entity: record.entity.0,
            polymer_kind: record.polymer_kind.native(),
        };
        self.0
            .push(native, checked_range(first_residue, end_residue)?)
            .map(PyChainIndex)
            .map_err(table_error)
    }

    fn residues(&self, chain: PyChainIndex) -> Option<[u32; 2]> {
        self.0.residues(chain.0).map(range_pair)
    }

    fn containing(&self, residue: u32) -> Option<PyChainIndex> {
        self.0.containing(residue).map(PyChainIndex)
    }

    fn label_asym_id(&self, chain: PyChainIndex) -> Option<PySymbolId> {
        self.0.label_asym_id(chain.0).map(PySymbolId)
    }

    fn auth_asym_id(&self, chain: PyChainIndex) -> Option<PySymbolId> {
        self.0.auth_asym_id(chain.0).map(PySymbolId)
    }

    fn entity(&self, chain: PyChainIndex) -> Option<PyEntityIndex> {
        self.0.entity(chain.0).map(PyEntityIndex)
    }

    fn polymer_kind(&self, chain: PyChainIndex) -> Option<PyPolymerKind> {
        self.0.polymer_kind(chain.0).map(Into::into)
    }

    fn rename(&mut self, chain: PyChainIndex, label: PySymbolId) -> PyResult<()> {
        self.0.rename(chain.0, label.0).map_err(table_error)
    }

    fn set_polymer_kind(&mut self, chain: PyChainIndex, kind: PyPolymerKind) -> PyResult<()> {
        self.0
            .set_polymer_kind(chain.0, kind.native())
            .map_err(table_error)
    }

    fn instances_of(&self, entity: PyEntityIndex) -> Vec<PyChainIndex> {
        self.0.instances_of(entity.0).map(PyChainIndex).collect()
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> Vec<PyChainIndex> {
        self.0.iter().map(PyChainIndex).collect()
    }
}

#[pymethods]
impl PyResidueRecord {
    #[new]
    #[pyo3(signature = (label_comp_id, auth_comp_id=None, label_seq_id=None, auth_seq_id=None, ins_code=None, het=false))]
    fn new(
        label_comp_id: PySymbolId,
        auth_comp_id: Option<PySymbolId>,
        label_seq_id: Option<i32>,
        auth_seq_id: Option<i32>,
        ins_code: Option<PySymbolId>,
        het: bool,
    ) -> Self {
        Self {
            label_comp_id,
            auth_comp_id: PyOptionalSymbol::new(auth_comp_id),
            label_seq_id: PyOptionalI32::new(label_seq_id),
            auth_seq_id: PyOptionalI32::new(auth_seq_id),
            ins_code: PyOptionalSymbol::new(ins_code),
            het,
        }
    }

    #[getter]
    fn label_comp_id(&self) -> PySymbolId {
        self.label_comp_id
    }

    #[getter]
    fn auth_comp_id(&self) -> PyOptionalSymbol {
        self.auth_comp_id
    }

    #[getter]
    fn label_seq_id(&self) -> PyOptionalI32 {
        self.label_seq_id
    }

    #[getter]
    fn auth_seq_id(&self) -> PyOptionalI32 {
        self.auth_seq_id
    }

    #[getter]
    fn ins_code(&self) -> PyOptionalSymbol {
        self.ins_code
    }

    #[getter]
    fn het(&self) -> bool {
        self.het
    }
}

#[pymethods]
impl PyResidueTable {
    #[new]
    fn new() -> Self {
        Self(pdbiox::core::topology::ResidueTable::default())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn push(
        &mut self,
        record: PyResidueRecord,
        first_atom: u32,
        end_atom: u32,
    ) -> PyResult<PyResidueIndex> {
        let native = pdbiox::core::topology::ResidueRecord {
            label_comp_id: record.label_comp_id.0,
            auth_comp_id: record.auth_comp_id.0,
            label_seq_id: record.label_seq_id.0,
            auth_seq_id: record.auth_seq_id.0,
            ins_code: record.ins_code.0,
            het: record.het,
        };
        self.0
            .push(native, checked_range(first_atom, end_atom)?)
            .map(PyResidueIndex)
            .map_err(table_error)
    }

    fn set_atoms(
        &mut self,
        residue: PyResidueIndex,
        first_atom: u32,
        end_atom: u32,
    ) -> PyResult<()> {
        self.0
            .set_atoms(residue.0, checked_range(first_atom, end_atom)?)
            .map_err(table_error)
    }

    fn atoms(&self, residue: PyResidueIndex) -> Option<[u32; 2]> {
        self.0.atoms(residue.0).map(range_pair)
    }

    fn containing(&self, atom: u32) -> Option<PyResidueIndex> {
        self.0.containing(atom).map(PyResidueIndex)
    }

    fn label_comp_id(&self, residue: PyResidueIndex) -> Option<PySymbolId> {
        self.0.label_comp_id(residue.0).map(PySymbolId)
    }

    fn auth_comp_id(&self, residue: PyResidueIndex) -> Option<PySymbolId> {
        self.0.auth_comp_id(residue.0).map(PySymbolId)
    }

    fn label_seq_id(&self, residue: PyResidueIndex) -> Option<i32> {
        self.0.label_seq_id(residue.0)
    }

    fn auth_seq_id(&self, residue: PyResidueIndex) -> Option<i32> {
        self.0.auth_seq_id(residue.0)
    }

    fn ins_code(&self, residue: PyResidueIndex) -> Option<PySymbolId> {
        self.0.ins_code(residue.0).map(PySymbolId)
    }

    fn is_het(&self, residue: PyResidueIndex) -> bool {
        self.0.is_het(residue.0)
    }
}

#[pymethods]
impl PyEntityTable {
    #[new]
    fn new() -> Self {
        Self(pdbiox::core::topology::EntityTable::default())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn push(
        &mut self,
        id: PyOptionalSymbol,
        kind: PyEntityKind,
        description: PyOptionalSymbol,
        canonical_sequence: Vec<PySymbolId>,
    ) -> PyResult<PyEntityIndex> {
        self.push_values(id, kind, description, canonical_sequence)
    }

    fn kind(&self, entity: PyEntityIndex) -> Option<PyEntityKind> {
        self.0.kind(entity.0).map(Into::into)
    }

    fn id(&self, entity: PyEntityIndex) -> Option<PySymbolId> {
        self.0.id(entity.0).map(PySymbolId)
    }

    fn find_by_id(&self, id: PySymbolId) -> Option<PyEntityIndex> {
        self.0.find_by_id(id.0).map(PyEntityIndex)
    }

    fn description(&self, entity: PyEntityIndex) -> Option<PySymbolId> {
        self.0.description(entity.0).map(PySymbolId)
    }

    fn canonical_sequence(&self, entity: PyEntityIndex) -> Vec<PySymbolId> {
        self.0
            .canonical_sequence(entity.0)
            .iter()
            .copied()
            .map(PySymbolId)
            .collect()
    }

    #[pyo3(name = "iter")]
    fn iter_values(&self) -> Vec<PyEntityIndex> {
        self.0.iter().map(PyEntityIndex).collect()
    }
}

impl PyEntityTable {
    pub(crate) fn push_values(
        &mut self,
        id: PyOptionalSymbol,
        kind: PyEntityKind,
        description: PyOptionalSymbol,
        canonical_sequence: Vec<PySymbolId>,
    ) -> PyResult<PyEntityIndex> {
        let sequence = canonical_sequence
            .into_iter()
            .map(|symbol| symbol.0)
            .collect::<Vec<_>>();
        let result = match id.0.get() {
            Some(id) => self.0.push(id, kind.native(), description.0, &sequence),
            None => self
                .0
                .push_without_id(kind.native(), description.0, &sequence),
        };
        result.map(PyEntityIndex).map_err(table_error)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyModelTable>()?;
    module.add_class::<PyChainRecord>()?;
    module.add_class::<PyChainTable>()?;
    module.add_class::<PyResidueRecord>()?;
    module.add_class::<PyResidueTable>()?;
    module.add_class::<PyEntityTable>()?;
    Ok(())
}
