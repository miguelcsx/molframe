//! CCD-topology perception for orbital atom types.

use super::{PeoeAtom, PeoeAtomType, PeoeBond, PeoeError};
use crate::Component;
use pdbiox_core::{BondOrder, Element};
use std::collections::BTreeMap;

pub(super) fn component_inputs(
    component: &Component,
) -> Result<(Vec<PeoeAtom>, Vec<PeoeBond>), PeoeError> {
    let mut names = BTreeMap::new();
    for (index, atom) in component.atoms.iter().enumerate() {
        if names.insert(atom.name.as_ref(), index).is_some() {
            return Err(PeoeError::DuplicateAtomName {
                atom_name: atom.name.clone(),
            });
        }
    }
    let mut adjacency = vec![Vec::new(); component.atoms.len()];
    let mut bonds = Vec::with_capacity(component.bonds.len());
    for bond in component.bonds.iter() {
        let atom_a = index_of(&names, &bond.atom_a)?;
        let atom_b = index_of(&names, &bond.atom_b)?;
        if atom_a == atom_b {
            return Err(PeoeError::SelfBond { atom: atom_a });
        }
        adjacency[atom_a].push((atom_b, bond.order));
        adjacency[atom_b].push((atom_a, bond.order));
        bonds.push(PeoeBond { atom_a, atom_b });
    }
    let mut atoms = Vec::with_capacity(component.atoms.len());
    for (index, atom) in component.atoms.iter().enumerate() {
        atoms.push(PeoeAtom {
            atom_type: perceive(index, component, &adjacency)?,
            formal_charge: f64::from(atom.charge),
        });
    }
    Ok((atoms, bonds))
}

fn index_of(names: &BTreeMap<&str, usize>, name: &str) -> Result<usize, PeoeError> {
    match names.get(name) {
        Some(index) => Ok(*index),
        None => Err(PeoeError::UnknownBondAtom {
            atom_name: name.into(),
        }),
    }
}

fn perceive(
    index: usize,
    component: &Component,
    adjacency: &[Vec<(usize, BondOrder)>],
) -> Result<PeoeAtomType, PeoeError> {
    let atom = &component.atoms[index];
    let aromatic = atom.aromatic
        || adjacency[index]
            .iter()
            .any(|(_, order)| *order == BondOrder::Aromatic);
    let highest = match adjacency[index]
        .iter()
        .map(|(_, order)| bond_rank(*order))
        .max()
    {
        Some(rank) => rank,
        None => 1,
    };
    let atom_type = match atom.element.atomic_number() {
        1 => PeoeAtomType::H,
        6 => orbital(
            highest,
            aromatic,
            PeoeAtomType::CSp3,
            PeoeAtomType::CSp2,
            PeoeAtomType::CSp,
        ),
        7 => orbital(
            highest,
            aromatic,
            PeoeAtomType::NSp3,
            PeoeAtomType::NSp2,
            PeoeAtomType::NSp,
        ),
        8 if highest >= 2 || aromatic => PeoeAtomType::OSp2,
        8 => PeoeAtomType::OSp3,
        9 => PeoeAtomType::FSp3,
        17 => PeoeAtomType::ClSp3,
        35 => PeoeAtomType::BrSp3,
        53 => PeoeAtomType::ISp3,
        16 => sulfur(index, component, adjacency, highest, aromatic),
        15 if highest >= 2 || aromatic => PeoeAtomType::PSp2,
        15 => PeoeAtomType::PSp3,
        14 => orbital(
            highest,
            aromatic,
            PeoeAtomType::SiSp3,
            PeoeAtomType::SiSp2,
            PeoeAtomType::SiSp,
        ),
        5 if highest >= 2 || aromatic => PeoeAtomType::BSp2,
        5 => PeoeAtomType::BSp3,
        4 if highest >= 2 || aromatic => PeoeAtomType::BeSp2,
        4 => PeoeAtomType::BeSp3,
        12 => orbital(
            highest,
            aromatic,
            PeoeAtomType::MgSp3,
            PeoeAtomType::MgSp2,
            PeoeAtomType::MgSp,
        ),
        13 if highest >= 2 || aromatic => PeoeAtomType::AlSp2,
        13 => PeoeAtomType::AlSp3,
        _ => {
            return Err(PeoeError::UnsupportedAtom {
                atom: index,
                element: atom.element,
                environment: "unparameterised element",
            });
        }
    };
    Ok(atom_type)
}

const fn orbital(
    rank: u8,
    aromatic: bool,
    sp3: PeoeAtomType,
    sp2: PeoeAtomType,
    sp: PeoeAtomType,
) -> PeoeAtomType {
    if rank >= 3 {
        sp
    } else if rank >= 2 || aromatic {
        sp2
    } else {
        sp3
    }
}

fn sulfur(
    index: usize,
    component: &Component,
    adjacency: &[Vec<(usize, BondOrder)>],
    highest: u8,
    aromatic: bool,
) -> PeoeAtomType {
    if aromatic
        || highest == 2
            && adjacency[index]
                .iter()
                .all(|(other, _)| component.atoms[*other].element != Element::OXYGEN)
    {
        return PeoeAtomType::SSp2;
    }
    let oxygen_count = adjacency[index]
        .iter()
        .filter(|(other, _)| component.atoms[*other].element == Element::OXYGEN)
        .count();
    match oxygen_count {
        1 => PeoeAtomType::SO,
        2.. => PeoeAtomType::SO2,
        0 => PeoeAtomType::SSp3,
    }
}

const fn bond_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Triple | BondOrder::Quadruple => 3,
        BondOrder::Double | BondOrder::Aromatic => 2,
        BondOrder::Single | BondOrder::Polymeric | BondOrder::Unknown => 1,
    }
}
