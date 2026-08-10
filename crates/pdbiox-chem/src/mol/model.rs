//! Lossless-enough `CTfile` records and writer failures.

use pdbiox_core::element::Element;

/// One atom of a small molecule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolAtom {
    /// Chemical element.
    pub element: Element,
    /// Cartesian position in angstrom.
    pub position: [f32; 3],
}

/// One zero-based molecular bond.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MolBond {
    /// First endpoint.
    pub first: usize,
    /// Second endpoint.
    pub second: usize,
    /// `CTfile` order: 1, 2, 3, or 4 for aromatic.
    pub order: u8,
}

/// Format-neutral molecular graph.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Molecule {
    /// Atoms in file order.
    pub atoms: Vec<MolAtom>,
    /// Bonds in file order.
    pub bonds: Vec<MolBond>,
}

/// MDL dialect used by a record.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MolVersion {
    /// Fixed-column V2000.
    #[default]
    V2000,
    /// Tagged V3000.
    V3000,
}

/// Optional `CTfile` attributes for one atom.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MolAtomMetadata {
    /// Formal charge, when declared.
    pub formal_charge: Option<i8>,
    /// Mass number, when declared.
    pub isotope: Option<u16>,
    /// MDL tetrahedral parity or V3000 CFG value.
    pub stereo_parity: Option<u8>,
}

/// Optional `CTfile` attributes for one bond.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MolBondMetadata {
    /// MDL bond stereo or V3000 CFG value.
    pub stereo: Option<u8>,
}

/// One SDF data field, preserving name and complete value.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SdfProperty {
    /// Text between angle brackets in the property header.
    pub name: Box<str>,
    /// Property value without its terminating blank line.
    pub value: Box<str>,
}

/// One complete MOL block or SDF record.
#[derive(Clone, Debug, PartialEq)]
pub struct MolRecord {
    /// First header line.
    pub name: Box<str>,
    /// Second header line.
    pub program: Box<str>,
    /// Third header line.
    pub comment: Box<str>,
    /// Parsed `CTfile` dialect.
    pub version: MolVersion,
    /// Molecular graph.
    pub molecule: Molecule,
    /// Atom attributes aligned with `molecule.atoms`.
    pub atom_metadata: Vec<MolAtomMetadata>,
    /// Bond attributes aligned with `molecule.bonds`.
    pub bond_metadata: Vec<MolBondMetadata>,
    /// SDF fields in source order.
    pub properties: Vec<SdfProperty>,
}

/// Malformed or unrepresentable MOL/SDF data.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolError {
    /// Required syntax or a numeric field is malformed.
    Malformed,
    /// Parallel metadata does not match the molecular graph.
    MetadataLength,
    /// A bond endpoint lies outside the atom table.
    BondEndpoint,
    /// A value cannot be represented by the selected dialect.
    Unrepresentable(MolVersion),
}

impl std::fmt::Display for MolError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Malformed => formatter.write_str("malformed MOL/SDF record"),
            Self::MetadataLength => {
                formatter.write_str("MOL/SDF metadata length does not match atoms or bonds")
            }
            Self::BondEndpoint => {
                formatter.write_str("MOL/SDF bond endpoint is outside the atom table")
            }
            Self::Unrepresentable(version) => {
                write!(formatter, "value cannot be represented by {version:?}")
            }
        }
    }
}

impl std::error::Error for MolError {}
