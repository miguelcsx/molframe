//! Final covalent connectivity refines isolated-component roles in O(bonds).

use super::atom_chemistry::AtomChemistry;
use molframe_core::{BondOrder, Element, structure::StructureData};

pub(super) fn carbonyl_pair(carbon: Element, other: Element, order: BondOrder) -> bool {
    carbon == Element::CARBON
        && matches!(other, Element::OXYGEN | Element::SULFUR)
        && order == BondOrder::Double
}

pub(super) fn refine_linked_acceptors(
    data: &StructureData,
    chemistry: &mut [Option<AtomChemistry>],
) {
    for bond in data.bonds.iter() {
        if bond.is_metal_coordination(&data.dictionary)
            || !data
                .atom(bond.atom_a)
                .zip(data.atom(bond.atom_b))
                .is_some_and(|(a, b)| super::alt_compatible(a, b))
        {
            continue;
        }
        for (carbon, nitrogen) in [(bond.atom_a, bond.atom_b), (bond.atom_b, bond.atom_a)] {
            if !chemistry
                .get(carbon.as_usize())
                .and_then(Option::as_ref)
                .is_some_and(|atom| atom.carbonyl)
            {
                continue;
            }
            if let Some(Some(atom)) = chemistry.get_mut(nitrogen.as_usize())
                && atom.element == Element::NITROGEN
            {
                atom.roles.acceptor = false;
            }
        }
    }
}
