//! Per-residue secondary-structure vocabulary shared by readers and analyses.

/// A residue-level structural assignment.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SecondaryStructure {
    /// The file and available analyses supplied no assignment.
    #[default]
    Unknown,
    /// A residue outside a helix or sheet.
    Coil,
    /// An alpha or other helical segment.
    Helix,
    /// A beta strand.
    Strand,
    /// A short hydrogen-bonded turn.
    Turn,
}
