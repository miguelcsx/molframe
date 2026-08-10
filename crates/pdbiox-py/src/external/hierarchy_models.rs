//! Hierarchical exports for Gemmi and `Bio.PDB`.

use super::common::{coordinate_array, project};
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::adapters::{ExportAtom, ExportResidue, TopologyExport};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const GEMMI_MODULE: &str = "gemmi";
const GEMMI_STRUCTURE: &str = "Structure";
const GEMMI_MODEL: &str = "Model";
const GEMMI_CHAIN: &str = "Chain";
const GEMMI_RESIDUE: &str = "Residue";
const GEMMI_ATOM: &str = "Atom";
const GEMMI_SEQUENCE_ID: &str = "SeqId";
const GEMMI_ELEMENT: &str = "Element";
const GEMMI_POSITION: &str = "Position";
const BIOPYTHON_BUILDER_MODULE: &str = "Bio.PDB.StructureBuilder";
const BIOPYTHON_BUILDER: &str = "StructureBuilder";
const BLANK_IDENTIFIER: &str = " ";
const HETEROGEN_FLAG: &str = "H";

#[pymethods]
impl PyStructure {
    /// Returns a Gemmi `Structure` for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_gemmi(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let structure_id = entry_id(self)?;
        let model_name = deposited_model_number(self, model)?.to_string();
        let gemmi = PyModule::import(py, GEMMI_MODULE)?;
        let structure = gemmi.getattr(GEMMI_STRUCTURE)?.call0()?;
        structure.setattr("name", structure_id)?;
        let model_object = gemmi.getattr(GEMMI_MODEL)?.call1((model_name,))?;
        populate_gemmi(&gemmi, &model_object, &export)?;
        structure.call_method1("add_model", (model_object,))?;
        Ok(structure.unbind())
    }

    /// Returns a `Bio.PDB` SMCRA structure for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_biopython(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let structure_id = entry_id(self)?;
        let module = PyModule::import(py, BIOPYTHON_BUILDER_MODULE)?;
        let builder = module.getattr(BIOPYTHON_BUILDER)?.call0()?;
        builder.call_method1("init_structure", (structure_id,))?;
        builder.call_method1("init_model", (model,))?;
        populate_biopython(py, &builder, &export)?;
        Ok(builder.call_method0("get_structure")?.unbind())
    }
}

fn entry_id(structure: &PyStructure) -> PyResult<&str> {
    structure
        .structure()
        .data()
        .entry
        .id
        .as_deref()
        .ok_or_else(|| PyValueError::new_err("structure has no deposited entry identifier"))
}

fn deposited_model_number(structure: &PyStructure, model: usize) -> PyResult<i32> {
    let ordinal = u32::try_from(model)
        .map_err(|_| PyValueError::new_err("model position exceeds the supported range"))?;
    structure
        .structure()
        .model(pdbiox::ModelIndex::new(ordinal))
        .and_then(pdbiox::ModelRef::number)
        .ok_or_else(|| PyValueError::new_err(format!("model {model} has no deposited number")))
}

fn populate_gemmi(
    gemmi: &Bound<'_, PyModule>,
    model: &Bound<'_, PyAny>,
    export: &TopologyExport,
) -> PyResult<()> {
    for chain in &export.chains {
        let chain_object = gemmi.getattr(GEMMI_CHAIN)?.call1((&chain.id,))?;
        for residue_index in chain.residues.clone() {
            let residue = &export.residues[residue_index];
            let residue_object = gemmi.getattr(GEMMI_RESIDUE)?.call0()?;
            residue_object.setattr("name", &residue.name)?;
            residue_object.setattr("seqid", gemmi_sequence_id(gemmi, residue)?)?;
            for atom_index in residue.atoms.clone() {
                residue_object
                    .call_method1("add_atom", (gemmi_atom(gemmi, &export.atoms[atom_index])?,))?;
            }
            chain_object.call_method1("add_residue", (residue_object,))?;
        }
        model.call_method1("add_chain", (chain_object,))?;
    }
    Ok(())
}

fn gemmi_sequence_id<'py>(
    gemmi: &Bound<'py, PyModule>,
    residue: &ExportResidue,
) -> PyResult<Bound<'py, PyAny>> {
    let number = residue_number(residue)?;
    match &residue.insertion_code {
        Some(code) => gemmi.getattr(GEMMI_SEQUENCE_ID)?.call1((number, code)),
        None => gemmi
            .getattr(GEMMI_SEQUENCE_ID)?
            .call1((number, BLANK_IDENTIFIER)),
    }
}

fn gemmi_atom<'py>(gemmi: &Bound<'py, PyModule>, atom: &ExportAtom) -> PyResult<Bound<'py, PyAny>> {
    let object = gemmi.getattr(GEMMI_ATOM)?.call0()?;
    object.setattr("name", &atom.name)?;
    object.setattr(
        "element",
        gemmi
            .getattr(GEMMI_ELEMENT)?
            .call1((&atom.element_symbol,))?,
    )?;
    let [x, y, z] = atom.position;
    object.setattr("pos", gemmi.getattr(GEMMI_POSITION)?.call1((x, y, z))?)?;
    object.setattr("occ", optional_float(atom.occupancy))?;
    object.setattr("b_iso", optional_float(atom.b_factor))?;
    if let Some(charge) = atom.formal_charge {
        object.setattr("charge", charge)?;
    }
    if let Some(serial) = atom.serial {
        object.setattr("serial", serial)?;
    }
    if let Some(altloc) = &atom.alternate_location {
        object.setattr("altloc", altloc)?;
    }
    Ok(object)
}

fn populate_biopython(
    py: Python<'_>,
    builder: &Bound<'_, PyAny>,
    export: &TopologyExport,
) -> PyResult<()> {
    let coordinates = coordinate_array(py, export)?;
    for chain in &export.chains {
        builder.call_method1("init_chain", (&chain.id,))?;
        builder.call_method1("init_seg", (BLANK_IDENTIFIER,))?;
        for residue_index in chain.residues.clone() {
            let residue = &export.residues[residue_index];
            builder.call_method1(
                "init_residue",
                (
                    &residue.name,
                    biopython_heterogen_flag(residue),
                    residue_number(residue)?,
                    blank_identifier(residue.insertion_code.as_deref()),
                ),
            )?;
            for atom_index in residue.atoms.clone() {
                let atom = &export.atoms[atom_index];
                builder.call_method1(
                    "init_atom",
                    (
                        &atom.name,
                        coordinates.get_item(atom_index)?,
                        atom.b_factor,
                        atom.occupancy,
                        blank_identifier(atom.alternate_location.as_deref()),
                        &atom.name,
                        atom.serial,
                        &atom.element_symbol,
                    ),
                )?;
            }
        }
    }
    Ok(())
}

fn residue_number(residue: &ExportResidue) -> PyResult<i32> {
    residue.number.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "residue in chain {} has no numeric identifier",
            residue.chain
        ))
    })
}

const fn biopython_heterogen_flag(residue: &ExportResidue) -> &'static str {
    if residue.is_heterogen {
        HETEROGEN_FLAG
    } else {
        BLANK_IDENTIFIER
    }
}

fn optional_float(value: Option<f32>) -> f32 {
    match value {
        Some(value) => value,
        None => f32::NAN,
    }
}

fn blank_identifier(value: Option<&str>) -> &str {
    match value {
        Some(value) => value,
        None => BLANK_IDENTIFIER,
    }
}
