//! `ParmEd` `Structure` export.

use super::common::{coordinate_array, project};
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::BondOrder;
use pdbiox::adapters::TopologyExport;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

const PARMED_MODULE: &str = "parmed";
const STRUCTURE_CLASS: &str = "Structure";
const ATOM_CLASS: &str = "Atom";
const BOND_CLASS: &str = "Bond";

#[pymethods]
impl PyStructure {
    /// Returns a `ParmEd` `Structure` for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_parmed(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let parmed = PyModule::import(py, PARMED_MODULE)?;
        let structure = parmed.getattr(STRUCTURE_CLASS)?.call0()?;
        let atom_type = parmed.getattr(ATOM_CLASS)?;
        let bond_type = parmed.getattr(BOND_CLASS)?;
        let atoms = add_atoms(py, &structure, &atom_type, &export)?;
        structure.setattr("coordinates", coordinate_array(py, &export)?)?;
        add_bonds(&structure, &bond_type, &export, &atoms)?;
        Ok(structure.unbind())
    }
}

fn add_atoms<'py>(
    py: Python<'py>,
    structure: &Bound<'py, PyAny>,
    atom_type: &Bound<'py, PyAny>,
    export: &TopologyExport,
) -> PyResult<Vec<Bound<'py, PyAny>>> {
    let mut objects = Vec::with_capacity(export.atoms.len());
    for atom in &export.atoms {
        let residue = &export.residues[atom.residue];
        let chain = &export.chains[residue.chain];
        let keywords = PyDict::new(py);
        keywords.set_item("name", &atom.name)?;
        keywords.set_item("atomic_number", atom.atomic_number)?;
        keywords.set_item("mass", atom.mass)?;
        if let Some(serial) = atom.serial {
            keywords.set_item("number", serial)?;
        }
        let object = atom_type.call((), Some(&keywords))?;
        set_optional_atom_fields(&object, atom)?;
        structure.call_method1(
            "add_atom",
            (
                &object,
                &residue.name,
                residue.number,
                &chain.id,
                residue.insertion_code.as_deref(),
            ),
        )?;
        objects.push(object);
    }
    Ok(objects)
}

fn set_optional_atom_fields(
    object: &Bound<'_, PyAny>,
    atom: &pdbiox::adapters::ExportAtom,
) -> PyResult<()> {
    if let Some(value) = atom.formal_charge {
        object.setattr("formal_charge", value)?;
    }
    if let Some(value) = atom.occupancy {
        object.setattr("occupancy", value)?;
    }
    if let Some(value) = atom.b_factor {
        object.setattr("bfactor", value)?;
    }
    if let Some(value) = &atom.alternate_location {
        object.setattr("altloc", value)?;
    }
    Ok(())
}

fn add_bonds<'py>(
    structure: &Bound<'py, PyAny>,
    bond_type: &Bound<'py, PyAny>,
    export: &TopologyExport,
    atoms: &[Bound<'py, PyAny>],
) -> PyResult<()> {
    let bonds = structure.getattr("bonds")?.cast_into::<PyList>()?;
    for bond in &export.bonds {
        let object = match numeric_order(bond.order) {
            Some(order) => bond_type.call1((&atoms[bond.atom_a], &atoms[bond.atom_b], order))?,
            None => bond_type.call1((&atoms[bond.atom_a], &atoms[bond.atom_b]))?,
        };
        bonds.append(object)?;
    }
    Ok(())
}

const fn numeric_order(order: BondOrder) -> Option<f32> {
    match order {
        BondOrder::Single => Some(1.0),
        BondOrder::Double => Some(2.0),
        BondOrder::Triple => Some(3.0),
        BondOrder::Quadruple => Some(4.0),
        BondOrder::Aromatic => Some(1.5),
        BondOrder::Polymeric | BondOrder::Unknown => None,
    }
}
