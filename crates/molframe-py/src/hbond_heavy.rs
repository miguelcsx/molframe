//! Hydrogen bonds for structures without hydrogens, from heavy-atom geometry.

use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use crate::table::{PyTable, TableBuilder};
use molframe::analysis::{
    HeavyAtomHydrogenBondOptions, HeavyRoleTable, heavy_atom_hydrogen_bonds_with_roles,
};
use pyo3::prelude::*;

/// Donor and acceptor roles come from residue and atom names, not from hydrogens.
/// Rows are sorted by `(donor, acceptor)`; `acceptor_angle` is NaN when the acceptor
/// has no covalent neighbour. Periodic images are not supported.
///
/// `role_table` (text) or `role_file` (a path) replaces the built-in typing table with one
/// the caller supplies; Rust parses and validates it. `Table.applied["roles"]` names the
/// table that applied: name, version, source and licence.
#[pyfunction]
#[pyo3(name = "heavy_atom_hydrogen_bonds", signature = (
    structure,
    *,
    max_distance=3.5,
    min_distance=2.2,
    min_acceptor_angle=Some(90.0),
    covalent_distance=1.9,
    include_same_residue=false,
    include_sulfur=false,
    role_table=None,
    role_file=None,
    backend="auto",
    context=None,
))]
#[allow(
    clippy::too_many_arguments,
    reason = "every option is an independent keyword of the Python signature"
)]
fn heavy_atom_hydrogen_bonds_py(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    min_distance: f32,
    min_acceptor_angle: Option<f64>,
    covalent_distance: f32,
    include_same_residue: bool,
    include_sulfur: bool,
    role_table: Option<&str>,
    role_file: Option<std::path::PathBuf>,
    backend: &str,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyTable> {
    let options = HeavyAtomHydrogenBondOptions {
        maximum_distance: max_distance,
        minimum_distance: min_distance,
        minimum_acceptor_angle_degrees: min_acceptor_angle,
        covalent_distance,
        include_same_residue_pairs: include_same_residue,
        include_sulfur,
        backend: crate::backend::parse(backend)?,
    };
    let supplied = match (role_table, role_file) {
        (Some(_), Some(_)) => {
            return Err(crate::error::failure(
                "give role_table or role_file, not both",
            ));
        }
        (Some(text), None) => Some(text.to_owned()),
        (None, Some(path)) => Some(std::fs::read_to_string(&path).map_err(|error| {
            crate::error::failure(format!("cannot read {}: {error}", path.display()))
        })?),
        (None, None) => None,
    };
    let custom = supplied
        .map(|text| HeavyRoleTable::parse(&text))
        .transpose()
        .map_err(crate::error::failure)?;
    let roles = custom.as_ref().unwrap_or_else(|| HeavyRoleTable::builtin());
    let identity = roles.identity().to_string();
    let source = structure.inner.clone();
    let bonds = crate::execution::run(py, context, |context| {
        heavy_atom_hydrogen_bonds_with_roles(source.engine(), roles, options, context)
    })?
    .map_err(crate::error::kernel)?;
    let indices = |column: &[molframe::AtomIndex]| -> Vec<u32> {
        column.iter().map(|index| index.get()).collect()
    };
    Ok(TableBuilder::new(py, bonds.len())
        .indices("donor", &indices(bonds.donor()))
        .indices("acceptor", &indices(bonds.acceptor()))
        .single("distance", bonds.distance())
        .double("acceptor_angle", bonds.acceptor_angle_degrees())
        .applied("roles", identity)
        .finish())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(heavy_atom_hydrogen_bonds_py, module)?)
}
