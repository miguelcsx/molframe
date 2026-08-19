//! `RDKit` molecular graph and conformer export.

use super::common::project;
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pdbiox::BondOrder;
use pdbiox::adapters::TopologyExport;
use pyo3::prelude::*;

const CHEM_MODULE: &str = "rdkit.Chem";
const GEOMETRY_MODULE: &str = "rdkit.Geometry";
const EDITABLE_MOLECULE_CLASS: &str = "RWMol";
const ATOM_CLASS: &str = "Atom";
const CONFORMER_CLASS: &str = "Conformer";
const POINT_CLASS: &str = "Point3D";
const BOND_TYPE_CLASS: &str = "BondType";

#[pymethods]
impl PyStructure {
    /// Returns an unsanitized `RDKit` molecule with one conformer.
    #[pyo3(signature = (namespace, *, model=0))]
    fn to_rdkit(
        &self,
        py: Python<'_>,
        namespace: PyNamespace,
        model: usize,
    ) -> PyResult<Py<PyAny>> {
        let export = project(self, model, namespace)?;
        let chemistry = PyModule::import(py, CHEM_MODULE)?;
        let geometry = PyModule::import(py, GEOMETRY_MODULE)?;
        let editable = chemistry.getattr(EDITABLE_MOLECULE_CLASS)?.call0()?;
        add_atoms(&editable, &chemistry, &export)?;
        add_bonds(&editable, &chemistry, &export)?;
        let molecule = editable.call_method0("GetMol")?;
        add_conformer(&molecule, &chemistry, &geometry, &export)?;
        Ok(molecule.unbind())
    }
}

fn add_atoms(
    editable: &Bound<'_, PyAny>,
    chemistry: &Bound<'_, PyModule>,
    export: &TopologyExport,
) -> PyResult<()> {
    let atom_type = chemistry.getattr(ATOM_CLASS)?;
    for atom in &export.atoms {
        let object = atom_type.call1((atom.atomic_number,))?;
        if let Some(charge) = atom.formal_charge {
            object.call_method1("SetFormalCharge", (charge,))?;
        }
        editable.call_method1("AddAtom", (object,))?;
    }
    Ok(())
}

fn add_bonds(
    editable: &Bound<'_, PyAny>,
    chemistry: &Bound<'_, PyModule>,
    export: &TopologyExport,
) -> PyResult<()> {
    let bond_types = chemistry.getattr(BOND_TYPE_CLASS)?;
    for bond in &export.bonds {
        let name = match bond.order {
            BondOrder::Single => "SINGLE",
            BondOrder::Double => "DOUBLE",
            BondOrder::Triple => "TRIPLE",
            BondOrder::Quadruple => "QUADRUPLE",
            BondOrder::Aromatic => "AROMATIC",
            BondOrder::Polymeric | BondOrder::Unknown => "UNSPECIFIED",
        };
        editable.call_method1(
            "AddBond",
            (bond.atom_a, bond.atom_b, bond_types.getattr(name)?),
        )?;
    }
    Ok(())
}

fn add_conformer(
    molecule: &Bound<'_, PyAny>,
    chemistry: &Bound<'_, PyModule>,
    geometry: &Bound<'_, PyModule>,
    export: &TopologyExport,
) -> PyResult<()> {
    let conformer = chemistry
        .getattr(CONFORMER_CLASS)?
        .call1((export.atoms.len(),))?;
    let point_type = geometry.getattr(POINT_CLASS)?;
    for (index, atom) in export.atoms.iter().enumerate() {
        let [x, y, z] = atom.position;
        let point = point_type.call1((x, y, z))?;
        conformer.call_method1("SetAtomPosition", (index, point))?;
    }
    molecule.call_method1("AddConformer", (conformer, true))?;
    Ok(())
}
