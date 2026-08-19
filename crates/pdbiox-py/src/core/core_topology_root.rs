//! Python handle for a complete hierarchy snapshot.

use crate::core_topology::{
    PyChainRecord, PyChainTable, PyEntityTable, checked_range, table_error,
};
use crate::core_topology::{PyModelTable, PyResidueRecord, PyResidueTable};
use crate::core_values::{PyOptionalSymbol, PySymbolId};
use crate::index::{PyChainIndex, PyEntityIndex, PyModelIndex, PyResidueIndex};
use crate::metadata::PyEntityKind;
use pyo3::prelude::*;

#[pyclass(name = "Topology", from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyTopology(pub(crate) pdbiox::Topology);

#[pymethods]
impl PyTopology {
    #[new]
    fn new() -> Self {
        Self(pdbiox::Topology::default())
    }

    fn atom_count(&self) -> u32 {
        self.0.atom_count()
    }

    #[getter]
    fn models(&self) -> PyModelTable {
        PyModelTable(self.0.models.clone())
    }

    #[getter]
    fn chains(&self) -> PyChainTable {
        PyChainTable(self.0.chains.clone())
    }

    #[getter]
    fn residues(&self) -> PyResidueTable {
        PyResidueTable(self.0.residues.clone())
    }

    #[getter]
    fn entities(&self) -> PyEntityTable {
        PyEntityTable(self.0.entities.clone())
    }

    fn push_model(
        &mut self,
        model_num: i32,
        first_chain: u32,
        end_chain: u32,
    ) -> PyResult<PyModelIndex> {
        self.0
            .models
            .push(model_num, checked_range(first_chain, end_chain)?)
            .map(PyModelIndex)
            .map_err(table_error)
    }

    fn push_chain(
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
            .chains
            .push(native, checked_range(first_residue, end_residue)?)
            .map(PyChainIndex)
            .map_err(table_error)
    }

    fn push_residue(
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
            .residues
            .push(native, checked_range(first_atom, end_atom)?)
            .map(PyResidueIndex)
            .map_err(table_error)
    }

    fn push_entity(
        &mut self,
        id: PyOptionalSymbol,
        kind: PyEntityKind,
        description: PyOptionalSymbol,
        canonical_sequence: Vec<PySymbolId>,
    ) -> PyResult<PyEntityIndex> {
        let mut entities = PyEntityTable(std::mem::take(&mut self.0.entities));
        let result = entities.push_values(id, kind, description, canonical_sequence);
        self.0.entities = entities.0;
        result
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTopology>()?;
    Ok(())
}
