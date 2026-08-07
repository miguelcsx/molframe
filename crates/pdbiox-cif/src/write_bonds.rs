//! Canonical `struct_conn` output.

use pdbiox_core::bond::{BondOrder, BondProvenance};
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

pub(super) fn write(out: &mut String, structure: &Structure) {
    if structure.data().bonds.is_empty() {
        return;
    }
    out.push_str(HEADER);
    for (position, bond) in structure.data().bonds.iter().enumerate() {
        let Some(atom_a) = structure.atom(bond.atom_a) else {
            continue;
        };
        let Some(atom_b) = structure.atom(bond.atom_b) else {
            continue;
        };
        let (Some(a), Some(b)) = (endpoint(structure, atom_a), endpoint(structure, atom_b)) else {
            continue;
        };
        let _ = writeln!(
            out,
            "bond{} covale {} {} {} {} {} {} {} {} {} {} {} {}",
            position + 1,
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
            provenance(bond.provenance),
        );
    }
    out.push_str("#\n");
}

struct Endpoint<'a> {
    chain: &'a str,
    seq: String,
    component: &'a str,
    atom: &'a str,
    alt: &'a str,
}

fn endpoint<'a>(structure: &'a Structure, atom: AtomRef<'a>) -> Option<Endpoint<'a>> {
    let residue = atom.residue()?;
    let chain = structure
        .data()
        .topology
        .chains
        .containing(residue.index().get())
        .and_then(|index| structure.chain(index))?;
    Some(Endpoint {
        chain: chain.label()?,
        seq: sequence(residue),
        component: atom.component_name()?,
        atom: atom.name()?,
        alt: match atom.alt_label() {
            Some(label) => label,
            None => ".",
        },
    })
}

fn sequence(residue: ResidueRef<'_>) -> String {
    match residue.label_seq_id() {
        Some(number) => number.to_string(),
        None => ".".to_owned(),
    }
}

fn order(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "SING",
        BondOrder::Double => "DOUB",
        BondOrder::Triple => "TRIP",
        BondOrder::Quadruple => "QUAD",
        BondOrder::Aromatic => "AROM",
        BondOrder::Polymeric | BondOrder::Unknown => "?",
    }
}

fn provenance(provenance: BondProvenance) -> &'static str {
    match provenance {
        BondProvenance::File => "pdbiox:file",
        BondProvenance::ChemicalComponentDictionary => "pdbiox:ccd",
        BondProvenance::InferredDistance => "pdbiox:inferred-distance",
        BondProvenance::User => "pdbiox:user",
    }
}
