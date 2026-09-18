//! Canonical `struct_conn` output with complete endpoint preflight.

use super::options::{CifWriteError, CifWriteOptions};
use super::value::{Quoted, quoted};
use molframe_core::bond::BondOrder;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use std::fmt::{self, Display, Formatter};

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
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<(), CifWriteError> {
    if structure.data().bonds.is_empty() {
        return Ok(());
    }
    let Some(connection_type) = options.connection_type_id() else {
        return Err(CifWriteError::MissingConnectionTypeId);
    };
    let connection_type = quoted(connection_type);
    let _ = out.write_str(HEADER);
    for (position, bond) in structure.data().bonds.iter().enumerate() {
        let atom_a = atom(structure, position, bond.atom_a.get())?;
        let atom_b = atom(structure, position, bond.atom_b.get())?;
        let a = endpoint(structure, atom_a, position, 1)?;
        let b = endpoint(structure, atom_b, position, 2)?;
        let _ = writeln!(
            out,
            "{} {} {} {} {} ?",
            position + 1,
            connection_type,
            a,
            b,
            order(bond.order),
        );
    }
    let _ = out.write_str("#\n");
    Ok(())
}

fn atom(structure: &Structure, bond: usize, atom: u32) -> Result<AtomRef<'_>, CifWriteError> {
    match structure.atom(molframe_core::index::AtomIndex::new(atom)) {
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
    let chain = chain_of(structure, residue).and_then(molframe_core::structure::ChainRef::label);
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

struct Endpoint<'a> {
    chain: Quoted<'a>,
    seq: Option<i32>,
    component: Quoted<'a>,
    atom: Quoted<'a>,
    alt: Option<Quoted<'a>>,
}

fn endpoint<'a>(
    structure: &'a Structure,
    atom: AtomRef<'a>,
    bond: usize,
    side: u8,
) -> Result<Endpoint<'a>, CifWriteError> {
    let Some(residue) = atom.residue() else {
        return Err(missing(bond, side, "residue"));
    };
    let Some(chain) =
        chain_of(structure, residue).and_then(molframe_core::structure::ChainRef::label)
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
        chain: quoted(chain),
        seq: residue.label_seq_id(),
        component: quoted(component),
        atom: quoted(name),
        alt: atom.alt_label().map(quoted),
    })
}

impl Display for Endpoint<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{} {} {} {} {}",
            self.chain,
            OptionalInteger(self.seq),
            self.component,
            self.atom,
            OptionalQuoted(self.alt.as_ref()),
        )
    }
}

struct OptionalInteger(Option<i32>);

impl Display for OptionalInteger {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(value) => Display::fmt(&value, formatter),
            None => formatter.write_str("."),
        }
    }
}

struct OptionalQuoted<'borrow, 'text>(Option<&'borrow Quoted<'text>>);

impl Display for OptionalQuoted<'_, '_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(value) => Display::fmt(value, formatter),
            None => formatter.write_str("."),
        }
    }
}

fn chain_of<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'_>,
) -> Option<molframe_core::structure::ChainRef<'a>> {
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
