//! Public inputs and errors for charge equalisation.

use pdbiox_core::Element;

/// Versioned parameter collection used by the equalisation kernel.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PeoeParameterProfile {
    /// Hybridisation-aware Gasteiger-Marsili parameters.
    #[default]
    GasteigerMarsili,
}

impl PeoeParameterProfile {
    /// Stable profile identifier recorded in provenance.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::GasteigerMarsili => "gasteiger-marsili",
        }
    }

    /// Version of the bundled numerical table and atom-typing rules.
    #[must_use]
    pub const fn version(self) -> &'static str {
        match self {
            Self::GasteigerMarsili => "1",
        }
    }
}

/// Parameterised orbital state of an atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeoeAtomType {
    /// Hydrogen.
    H,
    /// Tetrahedral carbon.
    CSp3,
    /// Trigonal carbon.
    CSp2,
    /// Linear carbon.
    CSp,
    /// Tetrahedral nitrogen.
    NSp3,
    /// Trigonal nitrogen.
    NSp2,
    /// Linear nitrogen.
    NSp,
    /// Tetrahedral oxygen.
    OSp3,
    /// Trigonal oxygen.
    OSp2,
    /// Fluorine.
    FSp3,
    /// Chlorine.
    ClSp3,
    /// Bromine.
    BrSp3,
    /// Iodine.
    ISp3,
    /// Tetrahedral sulfur.
    SSp3,
    /// Sulfur with one oxygen neighbour.
    SO,
    /// Sulfur with at least two oxygen neighbours.
    SO2,
    /// Trigonal sulfur.
    SSp2,
    /// Tetrahedral phosphorus.
    PSp3,
    /// Trigonal phosphorus.
    PSp2,
    /// Tetrahedral silicon.
    SiSp3,
    /// Trigonal silicon.
    SiSp2,
    /// Linear silicon.
    SiSp,
    /// Tetrahedral boron.
    BSp3,
    /// Trigonal boron.
    BSp2,
    /// Tetrahedral beryllium.
    BeSp3,
    /// Trigonal beryllium.
    BeSp2,
    /// Tetrahedral magnesium.
    MgSp3,
    /// Trigonal magnesium.
    MgSp2,
    /// Linear magnesium.
    MgSp,
    /// Tetrahedral aluminium.
    AlSp3,
    /// Trigonal aluminium.
    AlSp2,
}

/// One atom supplied to the allocation-free numerical kernel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeoeAtom {
    /// Perceived orbital state.
    pub atom_type: PeoeAtomType,
    /// Integral or caller-declared starting formal charge.
    pub formal_charge: f64,
}

/// One undirected bond supplied to the numerical kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PeoeBond {
    /// First atom index.
    pub atom_a: usize,
    /// Second atom index.
    pub atom_b: usize,
}

/// Explicit numerical definition of an equalisation run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeoeOptions {
    /// Number of charge-transfer iterations.
    pub passes: usize,
    /// Damping applied on the first iteration.
    pub initial_damping: f64,
    /// Multiplier applied to damping after every iteration.
    pub damping_factor: f64,
    /// Differences at or below this threshold transfer no charge.
    pub minimum_electronegativity_difference: f64,
    /// Versioned numerical parameter collection.
    pub profile: PeoeParameterProfile,
}

impl Default for PeoeOptions {
    fn default() -> Self {
        Self {
            passes: 12,
            initial_damping: 1.0,
            damping_factor: 0.5,
            minimum_electronegativity_difference: 1.0e-12,
            profile: PeoeParameterProfile::GasteigerMarsili,
        }
    }
}

/// A component cannot be charged without changing or guessing its chemistry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PeoeError {
    /// Iteration controls are non-finite or outside their valid range.
    InvalidOptions,
    /// Two component atoms have the same local identifier.
    DuplicateAtomName {
        /// Repeated component-local identifier.
        atom_name: Box<str>,
    },
    /// A component bond references an absent atom identifier.
    UnknownBondAtom {
        /// Missing component-local identifier.
        atom_name: Box<str>,
    },
    /// A bond connects an atom to itself.
    SelfBond {
        /// Self-linked atom index.
        atom: usize,
    },
    /// A low-level bond endpoint is outside the atom input.
    BondIndexOutsideAtomArray {
        /// Invalid endpoint.
        atom: usize,
        /// Number of supplied atoms.
        atom_count: usize,
    },
    /// No parameter exists for the perceived atomic environment.
    UnsupportedAtom {
        /// Atom index in the input.
        atom: usize,
        /// Chemical element.
        element: Element,
        /// Perception result that could not be parameterised.
        environment: &'static str,
    },
}

impl std::fmt::Display for PeoeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidOptions => formatter.write_str("invalid PEOE iteration options"),
            Self::DuplicateAtomName { atom_name } => {
                write!(formatter, "duplicate component atom name {atom_name}")
            }
            Self::UnknownBondAtom { atom_name } => {
                write!(
                    formatter,
                    "component bond references unknown atom {atom_name}"
                )
            }
            Self::SelfBond { atom } => write!(formatter, "PEOE bond {atom}-{atom} is self-linked"),
            Self::BondIndexOutsideAtomArray { atom, atom_count } => write!(
                formatter,
                "PEOE bond atom index {atom} is outside the {atom_count}-atom input"
            ),
            Self::UnsupportedAtom {
                atom,
                element,
                environment,
            } => write!(
                formatter,
                "atom {atom} ({element:?}) has unsupported PEOE environment {environment}"
            ),
        }
    }
}

impl std::error::Error for PeoeError {}
