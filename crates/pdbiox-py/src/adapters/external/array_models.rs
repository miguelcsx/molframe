//! Vectorised exports for array-oriented external object models.

use super::common::{coordinate_array, project};
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use numpy::IntoPyArray;
use numpy::ndarray::Array2;
use pdbiox::BondOrder;
use pdbiox::adapters::TopologyExport;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

const ASE_MODULE: &str = "ase";
const ASE_ATOMS_CLASS: &str = "Atoms";
const BIOTITE_MODULE: &str = "biotite.structure";
const BIOTITE_ARRAY_CLASS: &str = "AtomArray";
const BIOTITE_BOND_LIST_CLASS: &str = "BondList";

#[pymethods]
impl PyStructure {
    /// Returns a vectorised ASE `Atoms` object for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_ase(&self, py: Python<'_>, namespace: PyNamespace, model: usize) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let ase = PyModule::import(py, ASE_MODULE)?;
        let keywords = PyDict::new(py);
        keywords.set_item(
            "numbers",
            export
                .atoms
                .iter()
                .map(|atom| atom.atomic_number)
                .collect::<Vec<_>>(),
        )?;
        keywords.set_item("positions", coordinate_array(py, &export)?)?;
        keywords.set_item(
            "masses",
            export
                .atoms
                .iter()
                .map(|atom| atom.mass)
                .collect::<Vec<_>>(),
        )?;
        Ok(ase
            .getattr(ASE_ATOMS_CLASS)?
            .call((), Some(&keywords))?
            .unbind())
    }

    /// Returns a vectorised Biotite `AtomArray` for one model.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_biotite(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let biotite = PyModule::import(py, BIOTITE_MODULE)?;
        let array = biotite
            .getattr(BIOTITE_ARRAY_CLASS)?
            .call1((export.atoms.len(),))?;
        populate_biotite(py, &array, &export)?;
        let bonds = biotite
            .getattr(BIOTITE_BOND_LIST_CLASS)?
            .call1((export.atoms.len(), biotite_bonds(py, &export)?))?;
        array.setattr("bonds", bonds)?;
        Ok(array.unbind())
    }
}

fn populate_biotite(
    py: Python<'_>,
    array: &Bound<'_, PyAny>,
    export: &TopologyExport,
) -> PyResult<()> {
    array.setattr("coord", coordinate_array(py, export)?)?;
    array.setattr(
        "chain_id",
        atom_values(export, |_, _residue, chain| chain.id.clone()),
    )?;
    array.setattr("res_id", residue_numbers(export)?)?;
    array.setattr(
        "ins_code",
        atom_values(export, |_, residue, _| biotite_insertion_code(residue)),
    )?;
    array.setattr(
        "res_name",
        atom_values(export, |_, residue, _| residue.name.clone()),
    )?;
    array.setattr(
        "hetero",
        atom_values(export, |_, residue, _| residue.is_heterogen),
    )?;
    array.setattr(
        "atom_name",
        export
            .atoms
            .iter()
            .map(|atom| atom.name.clone())
            .collect::<Vec<_>>(),
    )?;
    array.setattr(
        "element",
        export
            .atoms
            .iter()
            .map(|atom| atom.element_symbol.clone())
            .collect::<Vec<_>>(),
    )?;
    Ok(())
}

fn biotite_insertion_code(residue: &pdbiox::adapters::ExportResidue) -> String {
    match &residue.insertion_code {
        Some(code) => code.clone(),
        None => String::new(),
    }
}

fn atom_values<T>(
    export: &TopologyExport,
    value: impl Fn(
        &pdbiox::adapters::ExportAtom,
        &pdbiox::adapters::ExportResidue,
        &pdbiox::adapters::ExportChain,
    ) -> T,
) -> Vec<T> {
    export
        .atoms
        .iter()
        .map(|atom| {
            let residue = &export.residues[atom.residue];
            value(atom, residue, &export.chains[residue.chain])
        })
        .collect()
}

fn residue_numbers(export: &TopologyExport) -> PyResult<Vec<i32>> {
    export
        .atoms
        .iter()
        .map(|atom| {
            export.residues[atom.residue].number.ok_or_else(|| {
                PyValueError::new_err(format!(
                    "residue {} has no numeric identifier",
                    atom.residue
                ))
            })
        })
        .collect()
}

fn biotite_bonds<'py>(
    py: Python<'py>,
    export: &TopologyExport,
) -> PyResult<Bound<'py, numpy::PyArray2<i64>>> {
    const BOND_COLUMNS: usize = 3;
    let mut values = Vec::with_capacity(export.bonds.len() * BOND_COLUMNS);
    for bond in &export.bonds {
        values.extend(biotite_bond(bond)?);
    }
    let array = Array2::from_shape_vec((export.bonds.len(), BOND_COLUMNS), values)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(array.into_pyarray(py))
}

fn biotite_bond(bond: &pdbiox::adapters::ExportBond) -> PyResult<[i64; 3]> {
    let kind = match bond.order {
        BondOrder::Unknown | BondOrder::Polymeric => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 9,
    };
    let atom_a = i64::try_from(bond.atom_a)
        .map_err(|_| PyValueError::new_err("bond endpoint exceeds the Biotite index range"))?;
    let atom_b = i64::try_from(bond.atom_b)
        .map_err(|_| PyValueError::new_err("bond endpoint exceeds the Biotite index range"))?;
    Ok([atom_a, atom_b, kind])
}
