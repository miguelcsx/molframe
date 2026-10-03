//! Python constructors for immutable molecular query expressions.

use crate::bindings::PyQuery;
use pyo3::prelude::*;

macro_rules! selector {
    ($rust:ident, $builder:ident) => {
        #[pyfunction]
        fn $rust() -> PyQuery {
            PyQuery::from_native(molframe::Query::from_builder(
                molframe::query::col::$builder(),
            ))
        }
    };
}

selector!(all, all);
selector!(none, none);
selector!(protein, is_protein);
selector!(backbone, backbone);
selector!(sidechain, sidechain);
selector!(nucleic, nucleic);
selector!(nucleic_backbone, nucleic_backbone);
selector!(nucleic_base, nucleic_base);
selector!(nucleic_sugar, nucleic_sugar);
selector!(water, water);
selector!(ions, ions);
selector!(lipids, lipids);
selector!(glycans, glycans);
selector!(hetero, hetero);
selector!(hydrogen, hydrogen);
selector!(polar_hydrogen, polar_hydrogen);
selector!(nonpolar_hydrogen, nonpolar_hydrogen);
selector!(heavy, heavy);
selector!(polymer, polymer);
selector!(ligands, ligands);
selector!(aromatic, aromatic);
selector!(helix, helix);
selector!(strand, strand);
selector!(sheet, sheet);
selector!(alpha_helix, alpha_helix);
selector!(helix_310, helix_310);
selector!(pi_helix, pi_helix);
selector!(polyproline, polyproline);
selector!(bridge, bridge);
selector!(turn, turn);
selector!(bend, bend);
selector!(coil, coil);

#[pyfunction]
fn chain(value: &str) -> PyQuery {
    PyQuery::from_native(molframe::Query::from_builder(
        molframe::query::col::chain().eq(value),
    ))
}

#[pyfunction]
fn residue(value: &str) -> PyQuery {
    PyQuery::from_native(molframe::Query::from_builder(
        molframe::query::col::resname().eq(value),
    ))
}

#[pyfunction]
fn atom(value: &str) -> PyQuery {
    PyQuery::from_native(molframe::Query::from_builder(
        molframe::query::col::name().eq(value),
    ))
}

#[pyfunction]
fn within(radius: f32, target: &PyQuery) -> PyQuery {
    PyQuery::from_native(molframe::Query::within(radius, target.native().clone()))
}

#[pyfunction]
fn residues_within(radius: f32, target: &PyQuery) -> PyQuery {
    PyQuery::from_native(molframe::Query::residues_within(
        radius,
        target.native().clone(),
    ))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    macro_rules! add {
        ($name:ident) => {
            module.add_function(wrap_pyfunction!($name, module)?)?;
        };
    }
    add!(all);
    add!(none);
    add!(protein);
    add!(backbone);
    add!(sidechain);
    add!(nucleic);
    add!(nucleic_backbone);
    add!(nucleic_base);
    add!(nucleic_sugar);
    add!(water);
    add!(ions);
    add!(lipids);
    add!(glycans);
    add!(hetero);
    add!(hydrogen);
    add!(polar_hydrogen);
    add!(nonpolar_hydrogen);
    add!(heavy);
    add!(polymer);
    add!(ligands);
    add!(aromatic);
    add!(helix);
    add!(strand);
    add!(sheet);
    add!(alpha_helix);
    add!(helix_310);
    add!(pi_helix);
    add!(polyproline);
    add!(bridge);
    add!(turn);
    add!(bend);
    add!(coil);
    add!(chain);
    add!(residue);
    add!(atom);
    add!(within);
    add!(residues_within);
    Ok(())
}

#[cfg(test)]
#[path = "selection_expr_tests.rs"]
mod tests;
