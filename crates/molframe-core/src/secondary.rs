//! Per-residue secondary-structure vocabulary shared by readers and analyses.

/// A residue-level structural assignment, in the DSSP vocabulary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SecondaryStructure {
    /// The file and available analyses supplied no assignment.
    #[default]
    Unknown,
    /// A residue outside every helix, ladder, turn and bend (DSSP blank).
    Coil,
    /// An α-helix, i → i+4 (DSSP `H`).
    AlphaHelix,
    /// A 3₁₀-helix, i → i+3 (DSSP `G`).
    ThreeTenHelix,
    /// A π-helix, i → i+5 (DSSP `I`).
    PiHelix,
    /// A helix whose class the source does not name as α, 3₁₀ or π.
    OtherHelix,
    /// An isolated β-bridge, a single bridge pair outside a ladder (DSSP `B`).
    BetaBridge,
    /// A residue in a β-ladder of two or more bridges (DSSP `E`).
    Strand,
    /// A hydrogen-bonded turn (DSSP `T`).
    Turn,
    /// A bend: Cα direction changes by more than 70° (DSSP `S`).
    Bend,
}

impl SecondaryStructure {
    /// Whether the state is any kind of helix.
    #[must_use]
    pub const fn is_helix(self) -> bool {
        matches!(
            self,
            Self::AlphaHelix | Self::ThreeTenHelix | Self::PiHelix | Self::OtherHelix
        )
    }

    /// The stable one-byte code exchanged with consumers such as the Python
    /// binding and viewers. Codes 0–4 predate the 10-state vocabulary and keep
    /// their meaning; the later states are appended.
    #[must_use]
    pub const fn code(self) -> u8 {
        match self {
            Self::Unknown => 0,
            Self::Coil => 1,
            Self::AlphaHelix => 2,
            Self::Strand => 3,
            Self::Turn => 4,
            Self::ThreeTenHelix => 5,
            Self::PiHelix => 6,
            Self::OtherHelix => 7,
            Self::BetaBridge => 8,
            Self::Bend => 9,
        }
    }
}

/// Where a residue's secondary-structure state came from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SecondarySource {
    /// Nothing assigned the state.
    #[default]
    None,
    /// The deposited file's helix and sheet records.
    File,
    /// Kabsch–Sander assignment from backbone hydrogen bonds.
    Dssp,
    /// The Cα-only fallback for a trace without backbone atoms.
    CaOnly,
}

impl SecondarySource {
    /// Precedence when two sources disagree; the higher rank wins.
    #[must_use]
    pub const fn rank(self) -> u8 {
        match self {
            Self::None => 0,
            Self::CaOnly => 1,
            Self::Dssp => 2,
            Self::File => 3,
        }
    }

    /// The stable one-byte code exchanged with consumers.
    #[must_use]
    pub const fn code(self) -> u8 {
        match self {
            Self::None => 0,
            Self::File => 1,
            Self::Dssp => 2,
            Self::CaOnly => 3,
        }
    }
}

/// One residue's state together with where it came from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct SecondaryAssignment {
    /// The assigned state.
    pub state: SecondaryStructure,
    /// The source that assigned it.
    pub source: SecondarySource,
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
