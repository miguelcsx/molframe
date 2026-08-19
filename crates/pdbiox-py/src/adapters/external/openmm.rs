//! `OpenMM` `Topology` and vectorised position export.

use super::common::{coordinate_array, project};
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::BondOrder;
use pdbiox::adapters::{ExportResidue, TopologyExport};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

const APP_MODULE: &str = "openmm.app";
const UNIT_MODULE: &str = "openmm.unit";
const TOPOLOGY_CLASS: &str = "Topology";
const ELEMENT_NAMESPACE: &str = "element";
const ELEMENT_CLASS: &str = "Element";
const ANGSTROM_UNIT: &str = "angstrom";

#[pymethods]
impl PyStructure {
    /// Returns `(openmm.app.Topology, positions_in_angstrom)` for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_openmm(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let app = PyModule::import(py, APP_MODULE)?;
        let topology = app.getattr(TOPOLOGY_CLASS)?.call0()?;
        let elements = app.getattr(ELEMENT_NAMESPACE)?.getattr(ELEMENT_CLASS)?;
        let chain_objects = add_chains(&topology, &export)?;
        let residue_objects = add_residues(py, &topology, &export, &chain_objects)?;
        let atom_objects = add_atoms(py, &topology, &elements, &export, &residue_objects)?;
        add_bonds(py, &topology, &export, &atom_objects)?;

        let unit = PyModule::import(py, UNIT_MODULE)?.getattr(ANGSTROM_UNIT)?;
        let positions = coordinate_array(py, &export)?.call_method1("__mul__", (unit,))?;
        Ok(PyTuple::new(py, [topology.as_any(), &positions])?
            .into_any()
            .unbind())
    }
}

fn add_chains<'py>(
    topology: &Bound<'py, PyAny>,
    export: &TopologyExport,
) -> PyResult<Vec<Bound<'py, PyAny>>> {
    export
        .chains
        .iter()
        .map(|chain| topology.call_method1("addChain", (&chain.id,)))
        .collect()
}

fn add_residues<'py>(
    py: Python<'py>,
    topology: &Bound<'py, PyAny>,
    export: &TopologyExport,
    chains: &[Bound<'py, PyAny>],
) -> PyResult<Vec<Bound<'py, PyAny>>> {
    export
        .residues
        .iter()
        .map(|residue| {
            let keywords = residue_keywords(py, residue)?;
            topology.call_method(
                "addResidue",
                (&residue.name, &chains[residue.chain]),
                Some(&keywords),
            )
        })
        .collect()
}

fn residue_keywords<'py>(py: Python<'py>, residue: &ExportResidue) -> PyResult<Bound<'py, PyDict>> {
    let keywords = PyDict::new(py);
    if let Some(number) = residue.number {
        keywords.set_item("id", number.to_string())?;
    }
    if let Some(insertion_code) = &residue.insertion_code {
        keywords.set_item("insertionCode", insertion_code)?;
    }
    Ok(keywords)
}

fn add_atoms<'py>(
    py: Python<'py>,
    topology: &Bound<'py, PyAny>,
    elements: &Bound<'py, PyAny>,
    export: &TopologyExport,
    residues: &[Bound<'py, PyAny>],
) -> PyResult<Vec<Bound<'py, PyAny>>> {
    export
        .atoms
        .iter()
        .map(|atom| {
            let element = elements.call_method1("getByAtomicNumber", (atom.atomic_number,))?;
            let keywords = PyDict::new(py);
            if let Some(serial) = atom.serial {
                keywords.set_item("id", serial.to_string())?;
            }
            topology.call_method(
                "addAtom",
                (&atom.name, element, &residues[atom.residue]),
                Some(&keywords),
            )
        })
        .collect()
}

fn add_bonds(
    py: Python<'_>,
    topology: &Bound<'_, PyAny>,
    export: &TopologyExport,
    atoms: &[Bound<'_, PyAny>],
) -> PyResult<()> {
    for bond in &export.bonds {
        let keywords = PyDict::new(py);
        if let Some(order) = numeric_order(bond.order) {
            keywords.set_item("order", order)?;
        }
        topology.call_method(
            "addBond",
            (&atoms[bond.atom_a], &atoms[bond.atom_b]),
            Some(&keywords),
        )?;
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
