//! Public topology records decoded from TPR.

/// Metadata carried by the TPR header.
#[derive(Clone, Debug, PartialEq)]
pub struct TprHeader {
    /// Producer version text.
    pub producer: String,
    /// Format version.
    pub format_version: i32,
    /// Format generation.
    pub generation: i32,
    /// Scalar precision in bytes.
    pub precision: u8,
    /// Atom count declared by the header.
    pub atom_count: usize,
    /// File tag.
    pub tag: String,
    /// Whether coordinates are present after the topology.
    pub has_coordinates: bool,
}

/// One residue in a decoded GROMACS topology.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TprResidue {
    /// Residue name.
    pub name: String,
    /// One-based residue number within the expanded system.
    pub number: i32,
    /// Parent molecule instance.
    pub molecule: usize,
    /// Parent molecule-type name.
    pub molecule_type: String,
}

/// One atom in a decoded GROMACS topology.
#[derive(Clone, Debug, PartialEq)]
pub struct TprAtom {
    /// Expanded zero-based atom position.
    pub index: usize,
    /// Atom name.
    pub name: String,
    /// Force-field atom type.
    pub atom_type: String,
    /// Parent residue position.
    pub residue: usize,
    /// Mass in unified atomic mass units.
    pub mass: f64,
    /// Partial charge in elementary-charge units.
    pub charge: f64,
    /// Atomic number, absent for virtual sites and coarse-grained particles.
    pub atomic_number: Option<u8>,
}

/// One expanded covalent or constrained bond.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TprBond {
    /// First endpoint.
    pub atom_a: usize,
    /// Second endpoint.
    pub atom_b: usize,
}

/// A complete expanded topology from one TPR file.
#[derive(Clone, Debug, PartialEq)]
pub struct TprTopology {
    /// Header metadata.
    pub header: TprHeader,
    /// Expanded residues.
    pub residues: Vec<TprResidue>,
    /// Expanded atoms.
    pub atoms: Vec<TprAtom>,
    /// Expanded bonds and constraints.
    pub bonds: Vec<TprBond>,
}
