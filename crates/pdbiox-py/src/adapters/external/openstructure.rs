//! `OpenStructure` entity export.

use super::common::project;
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::adapters::TopologyExport;
use pyo3::prelude::*;

const MOL_MODULE: &str = "ost.mol";
const GEOMETRY_MODULE: &str = "ost.geom";
const ENTITY_FACTORY: &str = "CreateEntity";
const EDIT_MODE_CLASS: &str = "EditMode";
const BUFFERED_EDIT: &str = "BUFFERED_EDIT";
const VECTOR_CLASS: &str = "Vec3";

#[pymethods]
impl PyStructure {
    /// Returns an `OpenStructure` molecular entity for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_openstructure(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let molecule = PyModule::import(py, MOL_MODULE)?;
        let geometry = PyModule::import(py, GEOMETRY_MODULE)?;
        let entity = molecule.getattr(ENTITY_FACTORY)?.call0()?;
        let edit_mode = molecule.getattr(EDIT_MODE_CLASS)?.getattr(BUFFERED_EDIT)?;
        let editor = entity.call_method1("EditXCS", (edit_mode,))?;
        populate_entity(&editor, &geometry, &export)?;
        editor.call_method0("UpdateICS")?;
        Ok(entity.unbind())
    }
}

fn populate_entity(
    editor: &Bound<'_, PyAny>,
    geometry: &Bound<'_, PyModule>,
    export: &TopologyExport,
) -> PyResult<()> {
    let mut atoms = Vec::with_capacity(export.atoms.len());
    for chain in &export.chains {
        let chain_object = editor.call_method1("InsertChain", (&chain.id,))?;
        for residue_index in chain.residues.clone() {
            let residue = &export.residues[residue_index];
            let residue_object = match residue.number {
                Some(number) => {
                    editor.call_method1("AppendResidue", (&chain_object, &residue.name, number))?
                }
                None => editor.call_method1("AppendResidue", (&chain_object, &residue.name))?,
            };
            for atom_index in residue.atoms.clone() {
                let atom = &export.atoms[atom_index];
                let [x, y, z] = atom.position;
                let vector = geometry.getattr(VECTOR_CLASS)?.call1((x, y, z))?;
                atoms.push(editor.call_method1(
                    "InsertAtom",
                    (
                        &residue_object,
                        &atom.name,
                        vector,
                        &atom.element_symbol,
                        optional_float(atom.occupancy),
                        optional_float(atom.b_factor),
                        residue.is_heterogen,
                    ),
                )?);
            }
        }
    }
    for bond in &export.bonds {
        editor.call_method1("Connect", (&atoms[bond.atom_a], &atoms[bond.atom_b]))?;
    }
    Ok(())
}

fn optional_float(value: Option<f32>) -> f32 {
    match value {
        Some(value) => value,
        None => f32::NAN,
    }
}
