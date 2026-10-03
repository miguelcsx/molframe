//! Existing glycosidic edges and conservative vacant-site spatial inference.

use super::rings::{compatible, hetero};
use super::types::{CarbohydrateLink, CarbohydrateReport};
use crate::grid::{CellGrid, cell_for};
use molframe_core::structure::AtomRef;
use molframe_core::{AtomIndex, BondAdjacency, BondProvenance, Element, Structure};

pub(super) fn collect_links(
    structure: &Structure,
    report: &mut CarbohydrateReport,
    fallback: bool,
) {
    let graph = structure.data().bonds.adjacency(structure.atom_count());
    let mut donors: Vec<_> = report
        .monosaccharides
        .iter()
        .enumerate()
        .filter_map(|(ring, sugar)| sugar.anomeric_atom.map(|atom| (atom, ring)))
        .collect();
    donors.sort_unstable();
    for bond in structure.data().bonds.iter() {
        let pair = donors
            .binary_search_by_key(&bond.atom_a, |(atom, _)| *atom)
            .ok()
            .map(|position| (donors[position].1, bond.atom_a, bond.atom_b))
            .or_else(|| {
                donors
                    .binary_search_by_key(&bond.atom_b, |(atom, _)| *atom)
                    .ok()
                    .map(|position| (donors[position].1, bond.atom_b, bond.atom_a))
            });
        let Some((donor, carbon, acceptor)) = pair else {
            continue;
        };
        let Some(atom) = structure.atom(acceptor) else {
            continue;
        };
        if atom
            .residue()
            .map(molframe_core::structure::ResidueRef::index)
            == Some(report.monosaccharides[donor].residue)
            || !atom.element().is_some_and(hetero)
            || !compatible(structure, carbon, acceptor)
        {
            continue;
        }
        // An existing distance-inferred edge remains inference, not a declaration.
        // Apply the same chemical/site and distance boundary to that provenance.
        if bond.provenance == BondProvenance::InferredDistance
            && (!vacant_donor(
                structure,
                graph,
                &report.monosaccharides[donor],
                carbon,
                Some(acceptor),
            ) || !eligible_acceptor(structure, graph, atom, report, Some(carbon))
                || !in_window(structure, carbon, acceptor))
        {
            continue;
        }
        push_link(
            structure,
            graph,
            report,
            donor,
            carbon,
            acceptor,
            bond.provenance,
        );
    }
    if fallback {
        infer_links(structure, graph, &donors, report);
    }
    report.links.sort_by_key(|link| {
        (
            link.donor,
            link.acceptor,
            link.donor_atom,
            link.acceptor_atom,
        )
    });
    report
        .terminal_links
        .sort_by_key(|link| (link.donor, link.acceptor_atom));
}

fn acceptor_ring(
    structure: &Structure,
    graph: &BondAdjacency,
    report: &CarbohydrateReport,
    atom: AtomIndex,
) -> Option<usize> {
    let residue = structure.atom(atom)?.residue()?.index();
    let start = report
        .monosaccharides
        .partition_point(|sugar| sugar.residue < residue);
    let end = report
        .monosaccharides
        .partition_point(|sugar| sugar.residue <= residue);
    report.monosaccharides[start..end]
        .iter()
        .enumerate()
        .find_map(|(local, sugar)| {
            let index = start + local;
            if sugar.residue != residue {
                return None;
            }
            if sugar.ring_atoms.contains(&atom) {
                return Some(index);
            }
            // O4 is one edge from the ring; O6 is two edges via the exocyclic C6.
            let connected = graph.neighbours(atom).iter().any(|next| {
                sugar.ring_atoms.contains(next)
                    || graph
                        .neighbours(*next)
                        .iter()
                        .any(|next| sugar.ring_atoms.contains(next))
            });
            connected.then_some(index)
        })
}

fn push_link(
    structure: &Structure,
    graph: &BondAdjacency,
    report: &mut CarbohydrateReport,
    donor: usize,
    donor_atom: AtomIndex,
    acceptor_atom: AtomIndex,
    provenance: BondProvenance,
) {
    let acceptor = acceptor_ring(structure, graph, report, acceptor_atom);
    let link = CarbohydrateLink {
        donor,
        acceptor,
        donor_atom,
        acceptor_atom,
        provenance,
    };
    if acceptor.is_some() {
        report.links.push(link);
    } else {
        report.terminal_links.push(link);
    }
}

#[derive(Clone, Copy)]
struct Site {
    atom: AtomIndex,
    donor: Option<usize>,
}

fn vacant_donor(
    structure: &Structure,
    graph: &BondAdjacency,
    sugar: &super::types::Monosaccharide,
    atom: AtomIndex,
    exclude: Option<AtomIndex>,
) -> bool {
    graph.neighbours(atom).iter().all(|other| {
        Some(*other) == exclude
            || sugar.ring_atoms.contains(other)
            || structure
                .atom(*other)
                .and_then(molframe_core::structure::AtomRef::element)
                .is_some_and(|element| !hetero(element))
    })
}

fn infer_links(
    structure: &Structure,
    graph: &BondAdjacency,
    donors: &[(AtomIndex, usize)],
    report: &mut CarbohydrateReport,
) {
    let mut sites = Vec::new();
    for (atom, donor) in donors {
        let sugar = &report.monosaccharides[*donor];
        if vacant_donor(structure, graph, sugar, *atom, None) {
            add_site(
                structure,
                &mut sites,
                Site {
                    atom: *atom,
                    donor: Some(*donor),
                },
            );
        }
    }
    if sites.is_empty() {
        return;
    }
    for atom in structure.data().atoms() {
        if eligible_acceptor(structure, graph, atom, report, None) {
            add_site(
                structure,
                &mut sites,
                Site {
                    atom: atom.index(),
                    donor: None,
                },
            );
        }
    }
    let grid = CellGrid::build(sites);
    let mut candidates = Vec::new();
    grid.for_each_cell(|own, neighbourhood| {
        for site in &grid.items()[own] {
            let Some(donor) = site.donor else {
                continue;
            };
            for range in neighbourhood {
                for acceptor in &grid.items()[range.clone()] {
                    if acceptor.donor.is_some() || !compatible(structure, site.atom, acceptor.atom)
                    {
                        continue;
                    }
                    let other_residue = structure
                        .atom(acceptor.atom)
                        .and_then(molframe_core::structure::AtomRef::residue)
                        .map(molframe_core::structure::ResidueRef::index);
                    if other_residue == Some(report.monosaccharides[donor].residue) {
                        continue;
                    }
                    if in_window(structure, site.atom, acceptor.atom) {
                        candidates.push((donor, site.atom, acceptor.atom));
                    }
                }
            }
        }
    });
    candidates.sort_unstable();
    candidates.dedup();
    // A distance cannot decide between competing glycosylation sites. Refuse
    // ambiguity rather than greedily committing the nearest plausible contact.
    let mut acceptors: Vec<_> = candidates.iter().map(|(_, _, atom)| *atom).collect();
    acceptors.sort_unstable();
    for (position, (donor, carbon, acceptor)) in candidates.iter().enumerate() {
        let unique_donor = (position == 0 || candidates[position - 1].0 != *donor)
            && (position + 1 == candidates.len() || candidates[position + 1].0 != *donor);
        let unique_acceptor = acceptors.partition_point(|atom| atom < acceptor) + 1
            == acceptors.partition_point(|atom| atom <= acceptor);
        if unique_donor && unique_acceptor {
            push_link(
                structure,
                graph,
                report,
                *donor,
                *carbon,
                *acceptor,
                BondProvenance::InferredDistance,
            );
        }
    }
}

fn add_site(structure: &Structure, sites: &mut Vec<(crate::grid::CellKey, Site)>, site: Site) {
    let Some(position) = structure
        .atom(site.atom)
        .and_then(molframe_core::structure::AtomRef::position)
    else {
        return;
    };
    if let Some(cell) = cell_for(position, 2.0) {
        sites.push((cell, site));
    }
}

fn eligible_acceptor(
    structure: &Structure,
    graph: &BondAdjacency,
    atom: AtomRef<'_>,
    report: &CarbohydrateReport,
    exclude: Option<AtomIndex>,
) -> bool {
    let Some(residue) = atom.residue() else {
        return false;
    };
    let heavy_count = graph
        .neighbours(atom.index())
        .iter()
        .filter(|index| {
            Some(**index) != exclude
                && structure
                    .atom(**index)
                    .and_then(molframe_core::structure::AtomRef::element)
                    .is_some_and(|element| !element.is_hydrogen())
        })
        .count();
    if heavy_count != 1 {
        return false;
    }
    if let Some(ring) = acceptor_ring(structure, graph, report, atom.index()) {
        return atom.element() == Some(Element::OXYGEN)
            && !report.monosaccharides[ring]
                .ring_atoms
                .contains(&atom.index());
    }
    matches!(
        (residue.name(), atom.name(), atom.element()),
        (Some("ASN"), Some("ND2"), Some(Element::NITROGEN))
            | (Some("SER"), Some("OG"), Some(Element::OXYGEN))
            | (Some("THR"), Some("OG1"), Some(Element::OXYGEN))
            | (Some("TYR"), Some("OH"), Some(Element::OXYGEN))
            | (Some("CYS"), Some("SG"), Some(Element::SULFUR))
    )
}

fn in_window(structure: &Structure, a: AtomIndex, b: AtomIndex) -> bool {
    let (Some(a), Some(b)) = (
        structure
            .atom(a)
            .and_then(molframe_core::structure::AtomRef::position),
        structure
            .atom(b)
            .and_then(molframe_core::structure::AtomRef::position),
    ) else {
        return false;
    };
    let distance: f32 = a.iter().zip(b).map(|(a, b)| (a - b) * (a - b)).sum();
    distance.is_finite() && (1.0..=4.0).contains(&distance)
}
