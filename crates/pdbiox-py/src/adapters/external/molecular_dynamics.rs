//! Vectorised `MDAnalysis` and `MDTraj` exports.

use super::common::{coordinate_array, project};
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::adapters::TopologyExport;
use pyo3::prelude::*;
use pyo3::types::PyDict;

const MDANALYSIS_MODULE: &str = "MDAnalysis";
const MDANALYSIS_UNIVERSE: &str = "Universe";
const MDTRAJ_MODULE: &str = "mdtraj";
const MDTRAJ_TOPOLOGY: &str = "Topology";
const MDTRAJ_TRAJECTORY: &str = "Trajectory";
const ANGSTROM_TO_NANOMETRE: f32 = 0.1;

#[pymethods]
impl PyStructure {
    /// Returns an in-memory `MDAnalysis` `Universe` for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_mdanalysis(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let module = PyModule::import(py, MDANALYSIS_MODULE)?;
        let keywords = PyDict::new(py);
        keywords.set_item("n_residues", export.residues.len())?;
        keywords.set_item("n_segments", export.chains.len())?;
        keywords.set_item(
            "atom_resindex",
            export
                .atoms
                .iter()
                .map(|atom| atom.residue)
                .collect::<Vec<_>>(),
        )?;
        keywords.set_item(
            "residue_segindex",
            export
                .residues
                .iter()
                .map(|residue| residue.chain)
                .collect::<Vec<_>>(),
        )?;
        keywords.set_item("trajectory", true)?;
        let universe = module.getattr(MDANALYSIS_UNIVERSE)?.call_method(
            "empty",
            (export.atoms.len(),),
            Some(&keywords),
        )?;
        populate_mdanalysis(&universe, &export)?;
        universe
            .getattr("atoms")?
            .setattr("positions", coordinate_array(py, &export)?)?;
        Ok(universe.unbind())
    }

    /// Returns an `MDTraj` `Trajectory` containing one frame.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_mdtraj(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let module = PyModule::import(py, MDTRAJ_MODULE)?;
        let topology = module.getattr(MDTRAJ_TOPOLOGY)?.call0()?;
        populate_mdtraj(&module, &topology, &export)?;
        let nanometres =
            coordinate_array(py, &export)?.call_method1("__mul__", (ANGSTROM_TO_NANOMETRE,))?;
        let frames = nanometres.call_method1("reshape", ((1, export.atoms.len(), 3),))?;
        Ok(module
            .getattr(MDTRAJ_TRAJECTORY)?
            .call1((frames, topology))?
            .unbind())
    }
}

fn populate_mdanalysis(universe: &Bound<'_, PyAny>, export: &TopologyExport) -> PyResult<()> {
    for (name, values) in [
        (
            "names",
            export
                .atoms
                .iter()
                .map(|atom| atom.name.clone())
                .collect::<Vec<_>>(),
        ),
        (
            "types",
            export
                .atoms
                .iter()
                .map(|atom| atom.element_symbol.clone())
                .collect::<Vec<_>>(),
        ),
        (
            "elements",
            export
                .atoms
                .iter()
                .map(|atom| atom.element_symbol.clone())
                .collect::<Vec<_>>(),
        ),
    ] {
        universe.call_method1("add_TopologyAttr", (name, values))?;
    }
    universe.call_method1(
        "add_TopologyAttr",
        (
            "masses",
            export
                .atoms
                .iter()
                .map(|atom| atom.mass)
                .collect::<Vec<_>>(),
        ),
    )?;
    universe.call_method1(
        "add_TopologyAttr",
        (
            "resnames",
            export
                .residues
                .iter()
                .map(|residue| residue.name.clone())
                .collect::<Vec<_>>(),
        ),
    )?;
    universe.call_method1(
        "add_TopologyAttr",
        ("resids", required_residue_numbers(export)?),
    )?;
    universe.call_method1(
        "add_TopologyAttr",
        (
            "segids",
            export
                .chains
                .iter()
                .map(|chain| chain.id.clone())
                .collect::<Vec<_>>(),
        ),
    )?;
    universe.call_method1("add_TopologyAttr", ("chainIDs", atom_chain_ids(export)))?;
    universe.call_method1(
        "add_TopologyAttr",
        (
            "bonds",
            export
                .bonds
                .iter()
                .map(|bond| (bond.atom_a, bond.atom_b))
                .collect::<Vec<_>>(),
        ),
    )?;
    Ok(())
}

fn populate_mdtraj(
    module: &Bound<'_, PyModule>,
    topology: &Bound<'_, PyAny>,
    export: &TopologyExport,
) -> PyResult<()> {
    let mut chains = Vec::with_capacity(export.chains.len());
    for _ in &export.chains {
        chains.push(topology.call_method0("add_chain")?);
    }
    let mut residues = Vec::with_capacity(export.residues.len());
    for residue in &export.residues {
        let keywords = PyDict::new(module.py());
        if let Some(number) = residue.number {
            keywords.set_item("resSeq", number)?;
        }
        residues.push(topology.call_method(
            "add_residue",
            (&residue.name, &chains[residue.chain]),
            Some(&keywords),
        )?);
    }
    let element_lookup = module
        .getattr("element")?
        .getattr("Element")?
        .getattr("getByAtomicNumber")?;
    let mut atoms = Vec::with_capacity(export.atoms.len());
    for atom in &export.atoms {
        let keywords = PyDict::new(module.py());
        if let Some(serial) = atom.serial {
            keywords.set_item("serial", serial)?;
        }
        atoms.push(topology.call_method(
            "add_atom",
            (
                &atom.name,
                element_lookup.call1((atom.atomic_number,))?,
                &residues[atom.residue],
            ),
            Some(&keywords),
        )?);
    }
    for bond in &export.bonds {
        topology.call_method1("add_bond", (&atoms[bond.atom_a], &atoms[bond.atom_b]))?;
    }
    Ok(())
}

fn required_residue_numbers(export: &TopologyExport) -> PyResult<Vec<i32>> {
    export
        .residues
        .iter()
        .enumerate()
        .map(|(index, residue)| {
            residue.number.ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "residue {index} has no numeric identifier"
                ))
            })
        })
        .collect()
}

fn atom_chain_ids(export: &TopologyExport) -> Vec<String> {
    export
        .atoms
        .iter()
        .map(|atom| {
            let residue = &export.residues[atom.residue];
            export.chains[residue.chain].id.clone()
        })
        .collect()
}
