//! Crystal contacts as placed copies of whole chains.
//!
//! A crystal neighbour search reports atoms. An analysis that runs over "the
//! structure and its crystal contacts" needs a structure, and a half residue is
//! not one: its missing atoms would read as modelling gaps. So contacts are
//! lifted to whole chains. Every chain of the deposited asymmetric unit that has
//! an atom within the radius of the unit, under some symmetry operation and
//! lattice translation, is placed once under that operation and translation, and
//! the unit itself is kept. The result is an [`AssemblyView`], so the same
//! materialisation, chain labelling and instance identifiers as a biological
//! assembly apply.

use crate::view::InstanceRecord;
use crate::{
    AssemblyView, CellTransform, CrystalNeighborOptions, SymmetryExt, SymmetryOperation,
    collect_crystal_neighbors,
};
use molframe_core::{
    ChainIndex, Code, Diagnostic, ExecutionContext, InstanceId, ModelIndex, Structure,
};
use molframe_geom::Rigid;
use std::collections::BTreeSet;
use std::ops::Range;

/// The deposited chains with the half-open atom range each covers.
#[must_use]
pub fn chain_atom_ranges(structure: &Structure) -> Vec<(ChainIndex, Range<u32>)> {
    structure
        .data()
        .topology
        .chains
        .iter()
        .map(|chain| (chain, crate::view::chain_atoms(structure, chain)))
        .collect()
}

/// The Cartesian rigid motion of a symmetry operation followed by a lattice
/// translation.
///
/// # Errors
///
/// Returns [`Code::E6015`] for an improper operation (a reflection or an
/// inversion): applying it to a chiral molecule would produce its mirror image,
/// which is not a crystal contact of that molecule.
pub fn operation_motion(
    transform: &CellTransform,
    operation: &SymmetryOperation,
    lattice: [i32; 3],
) -> Result<Rigid, Diagnostic> {
    let apply = |cartesian: [f64; 3]| -> [f64; 3] {
        let moved = operation.apply_fractional(transform.to_fractional(cartesian));
        let shifted = [
            moved[0] + f64::from(lattice[0]),
            moved[1] + f64::from(lattice[1]),
            moved[2] + f64::from(lattice[2]),
        ];
        transform.to_cartesian(shifted)
    };
    let origin = apply([0.0; 3]);
    let mut columns = [[0.0; 3]; 3];
    for (axis, column) in columns.iter_mut().enumerate() {
        let mut unit = [0.0; 3];
        unit[axis] = 1.0;
        let image = apply(unit);
        *column = [
            image[0] - origin[0],
            image[1] - origin[1],
            image[2] - origin[2],
        ];
    }
    // `columns[k]` is the image of the k-th basis vector: the transpose of the matrix.
    let rotation = [
        [columns[0][0], columns[1][0], columns[2][0]],
        [columns[0][1], columns[1][1], columns[2][1]],
        [columns[0][2], columns[1][2], columns[2][2]],
    ];
    let motion = Rigid::new(rotation, origin);
    let determinant = motion.determinant();
    if (determinant - 1.0).abs() > 1e-6 {
        return Err(Diagnostic::new(Code::E6015)
            .with_context("determinant", determinant.to_string())
            .with_context("operation", &*operation.id));
    }
    Ok(motion)
}

/// Whether two motions agree to well below any coordinate precision.
fn same_motion(left: &Rigid, right: &Rigid) -> bool {
    let close = |a: f64, b: f64| (a - b).abs() < 1e-6;
    left.rotation
        .iter()
        .flatten()
        .zip(right.rotation.iter().flatten())
        .all(|(a, b)| close(*a, *b))
        && left
            .translation
            .iter()
            .zip(&right.translation)
            .all(|(a, b)| close(*a, *b))
}

impl AssemblyView {
    /// A view that places whole chains under explicit rigid motions.
    ///
    /// Each pair is one instance, in the order given. This is the constructor
    /// for arrangements the file does not name, such as crystal contacts.
    ///
    /// # Errors
    ///
    /// Returns [`Code::E6013`] for a chain the structure does not have.
    pub fn from_placements(
        structure: &Structure,
        id: &str,
        placements: &[(ChainIndex, Rigid)],
    ) -> Result<Self, Diagnostic> {
        let ranges = chain_atom_ranges(structure);
        let mut transforms: Vec<Rigid> = Vec::new();
        let mut instances = Vec::with_capacity(placements.len());
        for (chain, motion) in placements {
            let atoms = ranges
                .iter()
                .find(|(candidate, _)| candidate == chain)
                .map(|(_, range)| range.clone())
                .ok_or_else(|| {
                    Diagnostic::new(Code::E6013).with_context("chain", chain.get().to_string())
                })?;
            let transform = if let Some(index) = transforms.iter().position(|known| known == motion)
            {
                index
            } else {
                transforms.push(*motion);
                transforms.len() - 1
            };
            let instance = u32::try_from(instances.len())
                .map_err(|_| Diagnostic::new(Code::E6011).with_context("assembly", id))?;
            instances.push(InstanceRecord {
                source_chain: *chain,
                transform,
                instance_id: InstanceId::new(instance),
                atoms,
            });
        }
        Ok(Self::assemble(structure, id, instances, transforms))
    }
}

/// The deposited unit plus every chain that touches it from a neighbouring cell,
/// on either side.
///
/// `radius` is the largest atom-to-atom distance, in ångström, at which a
/// symmetry mate counts as a contact. The first instances are the deposited
/// chains under the identity; the contacts follow in a fixed order (operation,
/// then lattice translation, then chain), so the result does not depend on the
/// order in which the neighbour search happened to report them.
///
/// # Errors
///
/// Returns diagnostics when the structure has no cell or no symmetry operators
/// ([`Code::E6016`]), when an operation is improper ([`Code::E6015`]), or when
/// the neighbour search fails or exceeds its limits.
pub fn crystal_contact_view(
    structure: &Structure,
    radius: f32,
    context: &ExecutionContext,
) -> Result<AssemblyView, Diagnostic> {
    // The search's own ceiling on candidate images is a resource limit, so it comes from
    // the execution context like every other.
    let symmetry = structure
        .symmetry_set()
        .filter(|set| !set.operations().is_empty())
        .ok_or_else(|| {
            Diagnostic::new(Code::E6016).with_context("reason", "no symmetry operators")
        })?;
    let cell = structure
        .data()
        .cell
        .ok_or_else(|| Diagnostic::new(Code::E5004))?;
    let transform = CellTransform::new(&cell)?;
    let neighbours = collect_crystal_neighbors(
        structure,
        symmetry,
        ModelIndex::new(0),
        f64::from(radius),
        CrystalNeighborOptions {
            candidate_limit: context.image_search_limit(),
            ..CrystalNeighborOptions::default()
        },
        context,
    )?;
    let ranges = chain_atom_ranges(structure);
    let chain_of = |atom: u32| -> Option<ChainIndex> {
        ranges
            .iter()
            .find(|(_, range)| range.contains(&atom))
            .map(|(chain, _)| *chain)
    };
    // The search reports each contact once, as a pair of atoms in the unit and in one
    // image. The environment of the unit also holds the image that touches it from the
    // other side: the inverse motion applied to the chain that was near. Both are kept.
    let mut contacts: BTreeSet<(usize, [i32; 3], u32, u32)> = BTreeSet::new();
    for neighbour in &neighbours {
        if !neighbour.is_symmetry_mate(symmetry) {
            continue;
        }
        if let (Some(image), Some(near)) = (
            chain_of(neighbour.image_atom.get()),
            chain_of(neighbour.source_atom.get()),
        ) {
            contacts.insert((
                neighbour.operation,
                neighbour.lattice,
                image.get(),
                near.get(),
            ));
        }
    }
    let mut placements: Vec<(ChainIndex, Rigid)> = ranges
        .iter()
        .map(|(chain, _)| (*chain, Rigid::IDENTITY))
        .collect();
    let mut place = |chain: u32, motion: Rigid| {
        let chain = ChainIndex::new(chain);
        if !placements
            .iter()
            .any(|(known, other)| *known == chain && same_motion(other, &motion))
        {
            placements.push((chain, motion));
        }
    };
    for (operation, lattice, image, near) in contacts {
        let Some(operation) = symmetry.operations().get(operation) else {
            continue;
        };
        let motion = operation_motion(&transform, operation, lattice)?;
        place(image, motion);
        place(near, motion.inverse());
    }
    AssemblyView::from_placements(structure, &format!("crystal-{radius}"), &placements)
}

#[cfg(test)]
#[path = "crystal_system_tests.rs"]
mod tests;
