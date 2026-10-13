//! Hydrogen bonds between heavy atoms, for structures that have no hydrogens.
//!
//! Crystal and many cryo-EM models stop at the heavy atoms, so the
//! hydrogen-based detector has nothing to measure. This mode types donors and
//! acceptors from residue and atom names (see `heavy_role_table`), and accepts a
//! donor–acceptor pair when
//!
//! - the heavy-atom distance lies in `[minimum_distance, maximum_distance]`
//!   (3.5 Å by default; 2.2 Å below it is a clash or a bond, not a contact);
//! - the atoms are not covalent neighbours or neighbours-of-neighbours (1-2 and
//!   1-3), which also removes `N(i+1)…O(i)` across a peptide bond, and, unless
//!   asked for, do not share a residue;
//! - the angle donor–acceptor–antecedent is at least 90° for every covalent
//!   neighbour of the acceptor (`McDonald` and Thornton's acceptor test).
//!
//! No hydrogen position is invented: records carry the two heavy atoms, their
//! distance and the smallest acceptor-side angle. Covalent neighbours come from
//! the bond table when present and from a 1.9 Å heavy-atom distance query always,
//! so the mode works on a bare PDB file.
//!
//! The pair search streams through the shared spatial reduction, so memory is
//! `O(atoms + accepted bonds)` and time `O(atoms · local density)`; output is
//! sorted by `(donor, acceptor)` and does not depend on thread count. Periodic
//! images are not supported: the options carry no periodic switch.

use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, Structure};
use molframe_core::{AtomIndex, ExecutionContext};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

use super::heavy_role_table::HeavyRoleTable;
use super::heavy_roles::assign_roles;
use super::{HydrogenBondError, altlocs_compatible, orientations};
use crate::numeric::f64_to_f32;

/// Geometric and typing policy for heavy-atom hydrogen bonds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeavyAtomHydrogenBondOptions {
    /// Largest donor–acceptor distance in ångström. Default 3.5.
    pub maximum_distance: f32,
    /// Smallest donor–acceptor distance in ångström. Default 2.2.
    pub minimum_distance: f32,
    /// Smallest donor–acceptor–antecedent angle in degrees, applied to every
    /// covalent neighbour of the acceptor. `None` records the angle without
    /// filtering on it. Default `Some(90.0)`.
    pub minimum_acceptor_angle_degrees: Option<f64>,
    /// Heavy-atom distance in ångström under which two atoms count as
    /// covalently bonded. Default 1.9.
    pub covalent_distance: f32,
    /// Keep pairs inside one residue. Default `false`.
    pub include_same_residue_pairs: bool,
    /// Type Cys `SG` (donor and acceptor) and Met `SD` (acceptor). Default
    /// `false`: sulfur hydrogen bonds are weak and ambiguous.
    pub include_sulfur: bool,
    /// Spatial implementation. Default [`SpatialBackend::Auto`].
    pub backend: SpatialBackend,
}

impl Default for HeavyAtomHydrogenBondOptions {
    fn default() -> Self {
        Self {
            maximum_distance: 3.5,
            minimum_distance: 2.2,
            minimum_acceptor_angle_degrees: Some(90.0),
            covalent_distance: 1.9,
            include_same_residue_pairs: false,
            include_sulfur: false,
            backend: SpatialBackend::Auto,
        }
    }
}

/// One heavy-atom hydrogen bond. No hydrogen is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeavyAtomHydrogenBond {
    /// Donor heavy atom.
    pub donor: AtomIndex,
    /// Acceptor heavy atom.
    pub acceptor: AtomIndex,
    /// Donor–acceptor distance in ångström.
    pub distance: f32,
    /// Smallest donor–acceptor–antecedent angle in degrees over the acceptor's
    /// covalent neighbours; `NaN` when the acceptor has none (a lone water).
    pub acceptor_angle_degrees: f64,
}

define_soa_table! {
    /// Native columnar storage for heavy-atom hydrogen bonds.
    pub struct HeavyAtomHydrogenBondTable for HeavyAtomHydrogenBond {
        /// Donor atom indices.
        donor: AtomIndex,
        /// Acceptor atom indices.
        acceptor: AtomIndex,
        /// Donor-to-acceptor distances.
        distance: f32,
        /// Smallest donor-acceptor-antecedent angles; `NaN` without antecedents.
        acceptor_angle_degrees: f64,
    }
}

/// Finds hydrogen bonds from heavy-atom geometry alone.
///
/// Hydrogens present in the structure are ignored. Atoms from incompatible
/// alternate locations are never paired, and a bonded neighbour in another
/// conformer is not an antecedent. Donor and acceptor typing comes from
/// `roles`, a table the caller selected by name; nothing is built in.
///
/// # Errors
///
/// [`HydrogenBondError::InvalidOptions`] for a non-finite or inconsistent
/// distance or angle, or a spatial-search failure.
pub fn heavy_atom_hydrogen_bonds(
    structure: &Structure,
    roles: &HeavyRoleTable,
    options: HeavyAtomHydrogenBondOptions,
    context: &ExecutionContext,
) -> Result<HeavyAtomHydrogenBondTable, HydrogenBondError> {
    validate(options)?;
    let roles = assign_roles(structure, roles, options.include_sulfur);
    let covalent = covalent_graph(structure, options, context)?;
    let positions = structure.positions();
    let residue_of = |atom: u32| {
        structure
            .data()
            .atom(AtomIndex::new(atom))
            .and_then(AtomRef::residue)
            .map(molframe_core::structure::ResidueRef::index)
    };
    let query = PairQuery {
        positions,
        left: &roles.donors,
        right: &roles.acceptors,
        cutoff: options.maximum_distance,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: None,
        context,
    };
    let minimum_squared = options.minimum_distance * options.minimum_distance;
    let parts = reduce_pairs_within_unsorted(
        &query,
        Vec::new,
        |out: &mut Vec<HeavyAtomHydrogenBond>, pair| {
            if pair.distance_squared < minimum_squared {
                return;
            }
            for (donor, acceptor) in
                orientations(pair.first, pair.second, &roles.donors, &roles.acceptors)
            {
                if donor == acceptor
                    || !altlocs_compatible(structure, donor, acceptor)
                    || covalent.within_two_bonds(donor, acceptor)
                    || (!options.include_same_residue_pairs
                        && residue_of(donor).is_some_and(|r| Some(r) == residue_of(acceptor)))
                {
                    continue;
                }
                let angle = acceptor_angle(positions, &covalent, donor, acceptor);
                if let (Some(limit), Some(angle)) = (options.minimum_acceptor_angle_degrees, angle)
                    && angle < limit
                {
                    continue;
                }
                out.push(HeavyAtomHydrogenBond {
                    donor: AtomIndex::new(donor),
                    acceptor: AtomIndex::new(acceptor),
                    distance: f64_to_f32(f64::from(pair.distance_squared).sqrt()),
                    acceptor_angle_degrees: match angle {
                        Some(angle) => angle,
                        None => f64::NAN,
                    },
                });
            }
        },
    )?;
    let mut bonds: Vec<HeavyAtomHydrogenBond> = parts.into_iter().flatten().collect();
    bonds.sort_by_key(|bond| (bond.donor.get(), bond.acceptor.get()));
    bonds.dedup_by_key(|bond| (bond.donor.get(), bond.acceptor.get()));
    Ok(bonds.into_iter().collect())
}

fn validate(options: HeavyAtomHydrogenBondOptions) -> Result<(), HydrogenBondError> {
    let angle_ok = options
        .minimum_acceptor_angle_degrees
        .is_none_or(|angle| angle.is_finite() && (0.0..=180.0).contains(&angle));
    if options.maximum_distance.is_finite()
        && options.minimum_distance.is_finite()
        && options.minimum_distance >= 0.0
        && options.minimum_distance < options.maximum_distance
        && options.covalent_distance.is_finite()
        && options.covalent_distance > 0.0
        && angle_ok
    {
        Ok(())
    } else {
        Err(HydrogenBondError::InvalidOptions)
    }
}

/// Smallest donor–acceptor–antecedent angle in degrees, `None` when the
/// acceptor has no other covalent neighbour.
fn acceptor_angle(
    positions: &[[f32; 3]],
    covalent: &Covalent,
    donor: u32,
    acceptor: u32,
) -> Option<f64> {
    let donor_position = *positions.get(donor as usize)?;
    let acceptor_position = *positions.get(acceptor as usize)?;
    covalent
        .neighbours(acceptor)
        .iter()
        .filter(|&&antecedent| antecedent != donor)
        .filter_map(|&antecedent| {
            let antecedent_position = *positions.get(antecedent as usize)?;
            molframe_geom::angle(donor_position, acceptor_position, antecedent_position)
        })
        .map(f64::to_degrees)
        .min_by(f64::total_cmp)
}

/// Heavy-atom covalent adjacency in compressed rows; neighbours are sorted.
struct Covalent {
    offsets: Vec<usize>,
    neighbours: Vec<u32>,
}

impl Covalent {
    fn neighbours(&self, atom: u32) -> &[u32] {
        let atom = atom as usize;
        match (self.offsets.get(atom), self.offsets.get(atom + 1)) {
            (Some(&start), Some(&end)) => match self.neighbours.get(start..end) {
                Some(row) => row,
                None => &[],
            },
            _ => &[],
        }
    }

    /// Whether the atoms are bonded or share a bonded neighbour.
    fn within_two_bonds(&self, first: u32, second: u32) -> bool {
        let (a, b) = (self.neighbours(first), self.neighbours(second));
        if a.binary_search(&second).is_ok() {
            return true;
        }
        let (mut i, mut j) = (0, 0);
        while let (Some(x), Some(y)) = (a.get(i), b.get(j)) {
            match x.cmp(y) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => return true,
            }
        }
        false
    }
}

fn covalent_graph(
    structure: &Structure,
    options: HeavyAtomHydrogenBondOptions,
    context: &ExecutionContext,
) -> Result<Covalent, HydrogenBondError> {
    let count = structure.atom_count();
    let is_heavy = |atom: u32| {
        !structure
            .data()
            .atom(AtomIndex::new(atom))
            .and_then(AtomRef::element)
            .is_some_and(molframe_core::Element::is_hydrogen)
    };
    let heavy = AtomSelection::from_sorted((0..count).filter(|&atom| is_heavy(atom)).collect());
    let query = PairQuery {
        positions: structure.positions(),
        left: &heavy,
        right: &heavy,
        cutoff: options.covalent_distance,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: None,
        context,
    };
    let parts =
        reduce_pairs_within_unsorted(&query, Vec::new, |out: &mut Vec<(u32, u32)>, pair| {
            if pair.first != pair.second && altlocs_compatible(structure, pair.first, pair.second) {
                out.push((pair.first, pair.second));
            }
        })?;
    let mut edges: Vec<(u32, u32)> = parts.into_iter().flatten().collect();
    let bonds = &structure.data().bonds;
    if bonds.is_available() {
        let adjacency = bonds.adjacency(count);
        for first in 0..count {
            for neighbour in adjacency.neighbours(AtomIndex::new(first)) {
                let second = neighbour.get();
                if first < second
                    && is_heavy(first)
                    && is_heavy(second)
                    && altlocs_compatible(structure, first, second)
                {
                    edges.push((first, second));
                }
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();
    Ok(compress(count as usize, &edges))
}

fn compress(count: usize, edges: &[(u32, u32)]) -> Covalent {
    let mut degree = vec![0usize; count + 1];
    for &(a, b) in edges {
        for atom in [a, b] {
            if let Some(slot) = degree.get_mut(atom as usize + 1) {
                *slot += 1;
            }
        }
    }
    for atom in 0..count {
        let previous = degree.get(atom).copied();
        if let (Some(previous), Some(slot)) = (previous, degree.get_mut(atom + 1)) {
            *slot += previous;
        }
    }
    let mut cursor = degree.clone();
    let mut neighbours = vec![0u32; edges.len() * 2];
    for &(a, b) in edges {
        for (from, to) in [(a, b), (b, a)] {
            if let Some(slot) = cursor.get_mut(from as usize)
                && let Some(cell) = neighbours.get_mut(*slot)
            {
                *cell = to;
                *slot += 1;
            }
        }
    }
    for atom in 0..count {
        if let (Some(&start), Some(&end)) = (degree.get(atom), degree.get(atom + 1))
            && let Some(row) = neighbours.get_mut(start..end)
        {
            row.sort_unstable();
        }
    }
    Covalent {
        offsets: degree,
        neighbours,
    }
}

#[cfg(test)]
#[path = "heavy_tests.rs"]
mod tests;
