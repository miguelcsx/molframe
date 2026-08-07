//! Internal-coordinate value types and rebuilding.

use crate::place_atom;
use pdbiox_core::{AtomIndex, Code, Diagnostic};
use pdbiox_geom::{angle, dihedral, distance};

/// Two bond lengths and their included angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hedron {
    /// Length `I–J` in ångström.
    pub first_length: f64,
    /// Angle `I–J–K` in radians.
    pub angle: f64,
    /// Length `J–K` in ångström.
    pub second_length: f64,
}

impl Hedron {
    /// Measures a hedron, or returns `None` for coincident arms.
    #[must_use]
    pub fn from_points(i: [f32; 3], j: [f32; 3], k: [f32; 3]) -> Option<Self> {
        Some(Self {
            first_length: distance(i, j),
            angle: angle(i, j, k)?,
            second_length: distance(j, k),
        })
    }
}

/// A four-atom internal coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dihedron {
    /// Hedron `I–J–K`.
    pub first: Hedron,
    /// Bond angle `J–K–L` in radians.
    pub angle: f64,
    /// Bond length `K–L` in ångström.
    pub length: f64,
    /// Torsion `I–J–K–L` in radians.
    pub torsion: f64,
}

impl Dihedron {
    /// Measures four Cartesian points, or returns `None` for degenerate frames.
    #[must_use]
    pub fn from_points(i: [f32; 3], j: [f32; 3], k: [f32; 3], l: [f32; 3]) -> Option<Self> {
        Some(Self {
            first: Hedron::from_points(i, j, k)?,
            angle: angle(j, k, l)?,
            length: distance(k, l),
            torsion: dihedral(i, j, k, l)?,
        })
    }
}

/// One atom placed from three previously available references.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalAtom {
    /// Atom reconstructed by this record.
    pub atom: AtomIndex,
    /// Reference atoms `I`, `J`, `K`.
    pub references: [AtomIndex; 3],
    /// Internal coordinate of `I–J–K–atom`.
    pub coordinate: Dihedron,
}

/// A deterministic internal-coordinate forest.
#[derive(Clone, Debug, PartialEq)]
pub struct InternalCoordinates {
    atom_count: usize,
    seeds: Box<[(AtomIndex, [f32; 3])]>,
    atoms: Box<[InternalAtom]>,
}

impl InternalCoordinates {
    pub(crate) fn new(
        atom_count: usize,
        seeds: Vec<(AtomIndex, [f32; 3])>,
        atoms: Vec<InternalAtom>,
    ) -> Self {
        Self {
            atom_count,
            seeds: seeds.into_boxed_slice(),
            atoms: atoms.into_boxed_slice(),
        }
    }

    /// Cartesian seed atoms required to orient disconnected or shallow branches.
    #[must_use]
    pub fn seeds(&self) -> &[(AtomIndex, [f32; 3])] {
        &self.seeds
    }

    /// Atoms represented purely by bond-angle-torsion coordinates.
    #[must_use]
    pub fn atoms(&self) -> &[InternalAtom] {
        &self.atoms
    }

    /// Rebuilds all recorded atoms; originally absent coordinates remain absent.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an invalid reference order or degenerate frame.
    pub fn rebuild(&self) -> Result<Vec<Option<[f32; 3]>>, Diagnostic> {
        let mut positions = vec![None; self.atom_count];
        for (atom, position) in &self.seeds {
            let slot = positions.get_mut(atom.as_usize()).ok_or_else(invariant)?;
            *slot = Some(*position);
        }
        for atom in &self.atoms {
            let [i, j, k] = atom.references.map(AtomIndex::as_usize);
            let references = [
                positions.get(i).copied().flatten(),
                positions.get(j).copied().flatten(),
                positions.get(k).copied().flatten(),
            ];
            let [Some(i), Some(j), Some(k)] = references else {
                return Err(invariant());
            };
            let coordinate = atom.coordinate;
            let position = place_atom(
                i,
                j,
                k,
                coordinate.length,
                coordinate.angle,
                coordinate.torsion,
            )
            .ok_or_else(|| Diagnostic::new(Code::E5002))?;
            let slot = positions
                .get_mut(atom.atom.as_usize())
                .ok_or_else(invariant)?;
            *slot = Some(position);
        }
        Ok(positions)
    }
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}
