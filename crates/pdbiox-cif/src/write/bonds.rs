//! Canonical `struct_conn` output with complete endpoint preflight.

use super::options::{CifWriteError, CifWriteOptions};
use super::value::quote_text;
use pdbiox_core::bond::BondOrder;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use std::fmt::Write as _;

const HEADER: &str = "loop_\n\
_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.pdbx_ptnr1_label_alt_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_ptnr2_label_alt_id\n\
_struct_conn.pdbx_value_order\n_struct_conn.details\n";

pub(super) fn preflight(
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<(), CifWriteError> {
    if structure.data().bonds.is_empty() {
        return Ok(());
    }
    if !options.generates_connection_ids() {
        return Err(CifWriteError::ConnectionIdsNotEnabled);
    }
    if options.connection_type_id().is_none_or(str::is_empty) {
        return Err(CifWriteError::MissingConnectionTypeId);
    }
    for (position, bond) in structure.data().bonds.iter().enumerate() {
        let atom_a = atom(structure, position, bond.atom_a.get())?;
        let atom_b = atom(structure, position, bond.atom_b.get())?;
        validate_endpoint(structure, atom_a, position, 1)?;
        validate_endpoint(structure, atom_b, position, 2)?;
    }
    Ok(())
}

pub(super) fn write(
    out: &mut String,
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<(), CifWriteError> {
    if structure.data().bonds.is_empty() {
        return Ok(());
    }
    let Some(connection_type) = options.connection_type_id() else {
        return Err(CifWriteError::MissingConnectionTypeId);
    };
    let connection_type = quote_text(connection_type);
    out.push_str(HEADER);
    for (position, bond) in structure.data().bonds.iter().enumerate() {
        let atom_a = atom(structure, position, bond.atom_a.get())?;
        let atom_b = atom(structure, position, bond.atom_b.get())?;
        let a = endpoint(structure, atom_a, position, 1)?;
        let b = endpoint(structure, atom_b, position, 2)?;
        let _ = writeln!(
            out,
            "{} {} {} {} {} {} {} {} {} {} {} {} {} ?",
            position + 1,
            connection_type,
            a.chain,
            a.seq,
            a.component,
            a.atom,
            a.alt,
            b.chain,
            b.seq,
            b.component,
            b.atom,
            b.alt,
            order(bond.order),
        );
    }
    out.push_str("#\n");
    Ok(())
}

fn atom(structure: &Structure, bond: usize, atom: u32) -> Result<AtomRef<'_>, CifWriteError> {
    match structure.atom(pdbiox_core::index::AtomIndex::new(atom)) {
        Some(value) => Ok(value),
        None => Err(CifWriteError::InvalidBondAtom { bond, atom }),
    }
}

fn validate_endpoint(
    structure: &Structure,
    atom: AtomRef<'_>,
    bond: usize,
    endpoint: u8,
) -> Result<(), CifWriteError> {
    let Some(residue) = atom.residue() else {
        return Err(missing(bond, endpoint, "residue"));
    };
    let chain = chain_of(structure, residue).and_then(pdbiox_core::structure::ChainRef::label);
    if chain.is_none_or(str::is_empty) {
        return Err(missing(bond, endpoint, "label_asym_id"));
    }
    if atom.component_name().is_none_or(str::is_empty) {
        return Err(missing(bond, endpoint, "label_comp_id"));
    }
    if atom.name().is_none_or(str::is_empty) {
        return Err(missing(bond, endpoint, "label_atom_id"));
    }
    Ok(())
}

const fn missing(bond: usize, endpoint: u8, field: &'static str) -> CifWriteError {
    CifWriteError::MissingBondField {
        bond,
        endpoint,
        field,
    }
}

struct Endpoint {
    chain: String,
    seq: String,
    component: String,
    atom: String,
    alt: String,
}

fn endpoint(
    structure: &Structure,
    atom: AtomRef<'_>,
    bond: usize,
    side: u8,
) -> Result<Endpoint, CifWriteError> {
    let Some(residue) = atom.residue() else {
        return Err(missing(bond, side, "residue"));
    };
    let Some(chain) =
        chain_of(structure, residue).and_then(pdbiox_core::structure::ChainRef::label)
    else {
        return Err(missing(bond, side, "label_asym_id"));
    };
    let Some(component) = atom.component_name() else {
        return Err(missing(bond, side, "label_comp_id"));
    };
    let Some(name) = atom.name() else {
        return Err(missing(bond, side, "label_atom_id"));
    };
    Ok(Endpoint {
        chain: quote_text(chain),
        seq: sequence(residue),
        component: quote_text(component),
        atom: quote_text(name),
        alt: match atom.alt_label() {
            Some(label) => quote_text(label),
            None => ".".to_owned(),
        },
    })
}

fn sequence(residue: ResidueRef<'_>) -> String {
    match residue.label_seq_id() {
        Some(number) => number.to_string(),
        None => ".".to_owned(),
    }
}

fn chain_of<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'_>,
) -> Option<pdbiox_core::structure::ChainRef<'a>> {
    let index = structure
        .data()
        .topology
        .chains
        .containing(residue.index().get())?;
    structure.chain(index)
}

const fn order(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "SING",
        BondOrder::Double => "DOUB",
        BondOrder::Triple => "TRIP",
        BondOrder::Quadruple => "QUAD",
        BondOrder::Aromatic => "AROM",
        BondOrder::Polymeric | BondOrder::Unknown => "?",
    }
}
