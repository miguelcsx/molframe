//! The repeating unit a scale fixture is built from, and where copies go.
//!
//! Scale is reached by tiling a real structure rather than by generating atoms
//! from a distribution. The tile carries real bond density, real residue
//! composition and real local neighbour statistics, because it is a real
//! structure; a point cloud drawn from a distribution carries none of them and
//! would exercise the coordinate path and nothing else.
//!
//! Two properties of the tiling are deliberate. Copies are placed on a cubic
//! lattice with a gap, so no copy overlaps another and the neighbour count per
//! atom stays what it is in the source. And each copy is rotated by its own
//! derived sequence, so the fixture does not present the same orientation to
//! every spatial query.
//!
//! What tiling does **not** reproduce is global heterogeneity: every copy has
//! identical chemistry, so a dictionary that grows with structure size is
//! under-tested and a per-component cache sees an unrealistically high hit
//! rate. Both flatter exactly the code paths a performance programme optimises,
//! which is why any number measured here is reported alongside the downloaded
//! real-structure targets rather than instead of them.

use super::seed::Seed;
use crate::numeric::{f64_to_f32, usize_to_f64};
use pdbiox_core::element::Element;
use pdbiox_core::structure::Structure;

/// Gap between adjacent lattice cells, as a fraction of the tile extent.
///
/// A tenth of the tile keeps copies clear of one another under any rotation
/// while leaving the overall density close to the source structure's.
const LATTICE_GAP_FRACTION: f64 = 0.1;

/// Default positional jitter, in ångström.
///
/// Below the three decimal places a deposited coordinate carries, so it
/// perturbs the spatial index without inventing chemistry the source does not
/// have.
const DEFAULT_JITTER: f64 = 0.05;

/// One atom of the repeating tile.
///
/// Identifier fields are positions in the tile's shared name table rather than
/// owned strings, because a tile has tens of thousands of atoms and a few
/// hundred distinct names.
#[derive(Clone, Copy, Debug)]
pub struct TileAtom {
    /// Position in the tile's name table of the atom name.
    pub atom_name: u32,
    /// Position in the tile's name table of the component identifier.
    pub comp_id: u32,
    /// Position in the tile's name table of the chain label.
    pub asym_id: u32,
    /// The residue's sequence position, where the source recorded one.
    pub seq_id: i32,
    /// The element.
    pub element: Element,
    /// The atom's position in the source structure.
    pub position: [f32; 3],
}

/// A rigid placement of one tile copy.
#[derive(Clone, Copy, Debug)]
pub struct Placement {
    rotation: [[f64; 3]; 3],
    translation: [f64; 3],
    jitter: f64,
}

impl Placement {
    /// Carries a tile position into this copy's frame.
    #[must_use]
    pub fn apply(&self, position: [f32; 3], seed: &mut Seed) -> [f32; 3] {
        let source = [
            f64::from(position[0]),
            f64::from(position[1]),
            f64::from(position[2]),
        ];
        let mut placed = [0.0_f64; 3];
        for (axis, slot) in placed.iter_mut().enumerate() {
            let row = self.rotation[axis];
            *slot = row[0].mul_add(
                source[0],
                row[1].mul_add(source[1], row[2].mul_add(source[2], self.translation[axis])),
            ) + seed.next_signed(self.jitter);
        }
        [
            f64_to_f32(placed[0]),
            f64_to_f32(placed[1]),
            f64_to_f32(placed[2]),
        ]
    }
}

/// A repeating unit and the lattice its copies are placed on.
#[derive(Clone, Debug)]
pub struct Tile {
    names: Vec<Box<str>>,
    atoms: Vec<TileAtom>,
    spacing: [f64; 3],
    seed: Seed,
    jitter: f64,
}

impl Tile {
    /// Extracts a tile from the first model of `structure`.
    #[must_use]
    pub fn from_structure(structure: &Structure, seed: Seed) -> Self {
        let mut names: Vec<Box<str>> = Vec::new();
        let mut atoms = Vec::new();
        let mut extent = [0.0_f32; 3];

        for chain in structure.data().chains() {
            let asym_id = intern(&mut names, chain.label());
            for residue in chain.residues() {
                let comp_id = intern(&mut names, residue.name());
                let seq_id = match residue.label_seq_id() {
                    Some(seq_id) => seq_id,
                    None => 1,
                };
                for atom in residue.atoms() {
                    let Some(position) = atom.position() else {
                        continue;
                    };
                    for (axis, slot) in extent.iter_mut().enumerate() {
                        *slot = slot.max(position[axis].abs());
                    }
                    atoms.push(TileAtom {
                        atom_name: intern(&mut names, atom.name()),
                        comp_id,
                        asym_id,
                        seq_id,
                        element: match atom.element() {
                            Some(element) => element,
                            None => Element::UNKNOWN,
                        },
                        position,
                    });
                }
            }
        }

        let spacing = spacing_from(extent);
        Self {
            names,
            atoms,
            spacing,
            seed,
            jitter: DEFAULT_JITTER,
        }
    }

    /// Sets the positional jitter applied to each placed atom, in ångström.
    #[must_use]
    pub const fn with_jitter(mut self, jitter: f64) -> Self {
        self.jitter = jitter;
        self
    }

    /// The atoms of one copy, in source order.
    #[must_use]
    pub fn atoms(&self) -> &[TileAtom] {
        &self.atoms
    }

    /// The tile's shared name table.
    #[must_use]
    pub fn names(&self) -> &[Box<str>] {
        &self.names
    }

    /// The text at `position` in the name table, or the empty string.
    #[must_use]
    pub fn name(&self, position: u32) -> &str {
        let Ok(index) = usize::try_from(position) else {
            return "";
        };
        match self.names.get(index) {
            Some(name) => name,
            None => "",
        }
    }

    /// How many copies are needed to reach at least `atoms` atoms.
    #[must_use]
    pub fn copies_for(&self, atoms: u64) -> u64 {
        let Ok(per_copy) = u64::try_from(self.atoms.len()) else {
            return 0;
        };
        if per_copy == 0 {
            return 0;
        }
        atoms.div_ceil(per_copy)
    }

    /// The number of atoms `copies` copies produce.
    #[must_use]
    pub fn atoms_in(&self, copies: u64) -> u64 {
        let Ok(per_copy) = u64::try_from(self.atoms.len()) else {
            return 0;
        };
        copies.saturating_mul(per_copy)
    }

    /// Where copy `ordinal` sits, given a total of `copies`.
    #[must_use]
    pub fn placement(&self, ordinal: u64, copies: u64) -> Placement {
        let mut seed = self.seed.derive(ordinal);
        let side = lattice_side(copies);
        let cell = lattice_cell(ordinal, side);
        let mut translation = [0.0_f64; 3];
        for (axis, slot) in translation.iter_mut().enumerate() {
            *slot = usize_to_f64(cell[axis]) * self.spacing[axis];
        }
        Placement {
            rotation: rotation_from(&mut seed),
            translation,
            jitter: self.jitter,
        }
    }

    /// A fresh derived sequence for the jitter of copy `ordinal`.
    #[must_use]
    pub fn jitter_seed(&self, ordinal: u64) -> Seed {
        self.seed.derive(ordinal ^ u64::MAX)
    }
}

/// Returns the position of `text` in `names`, appending it when it is new.
fn intern(names: &mut Vec<Box<str>>, text: Option<&str>) -> u32 {
    let text = match text {
        Some(text) => text,
        None => "",
    };
    if let Some(position) = names.iter().position(|name| name.as_ref() == text) {
        return narrow(position);
    }
    let position = narrow(names.len());
    names.push(text.into());
    position
}

/// Narrows a name-table position, saturating rather than wrapping.
///
/// A tile has a few hundred distinct names, so the saturating branch is
/// unreachable in practice and exists so the conversion needs no absence helper.
fn narrow(position: usize) -> u32 {
    match u32::try_from(position) {
        Ok(position) => position,
        Err(_) => u32::MAX,
    }
}

/// Lattice spacing that clears the tile's own extent under any rotation.
///
/// The bounding radius is used rather than the per-axis extent, because a
/// rotated copy can reach its longest diagonal along any axis.
fn spacing_from(extent: [f32; 3]) -> [f64; 3] {
    let radius = f64::from(extent[0])
        .hypot(f64::from(extent[1]))
        .hypot(f64::from(extent[2]));
    let step = radius
        .mul_add(2.0, radius * 2.0 * LATTICE_GAP_FRACTION)
        .max(1.0);
    [step; 3]
}

/// The cube side long enough to hold `copies` lattice sites.
fn lattice_side(copies: u64) -> u64 {
    let mut side = 1u64;
    while side.saturating_mul(side).saturating_mul(side) < copies.max(1) {
        side += 1;
    }
    side
}

/// The lattice cell copy `ordinal` occupies.
fn lattice_cell(ordinal: u64, side: u64) -> [usize; 3] {
    let side = side.max(1);
    let plane = side.saturating_mul(side);
    let x = ordinal / plane;
    let remainder = ordinal % plane;
    [widen(x), widen(remainder / side), widen(remainder % side)]
}

/// Widens a lattice coordinate, saturating rather than wrapping.
fn widen(coordinate: u64) -> usize {
    match usize::try_from(coordinate) {
        Ok(coordinate) => coordinate,
        Err(_) => usize::MAX,
    }
}

/// A rotation matrix from three angles drawn from `seed`.
fn rotation_from(seed: &mut Seed) -> [[f64; 3]; 3] {
    let alpha = seed.next_unit() * std::f64::consts::TAU;
    let beta = seed.next_unit() * std::f64::consts::TAU;
    let gamma = seed.next_unit() * std::f64::consts::TAU;

    let (sin_a, cos_a) = alpha.sin_cos();
    let (sin_b, cos_b) = beta.sin_cos();
    let (sin_g, cos_g) = gamma.sin_cos();

    [
        [
            cos_a * cos_b,
            cos_a.mul_add(sin_b * sin_g, -(sin_a * cos_g)),
            cos_a.mul_add(sin_b * cos_g, sin_a * sin_g),
        ],
        [
            sin_a * cos_b,
            sin_a.mul_add(sin_b * sin_g, cos_a * cos_g),
            sin_a.mul_add(sin_b * cos_g, -(cos_a * sin_g)),
        ],
        [-sin_b, cos_b * sin_g, cos_b * cos_g],
    ]
}

#[cfg(test)]
#[path = "tile_tests.rs"]
mod tests;
