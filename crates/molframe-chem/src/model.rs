//! Chemical component definitions independent of structure-local symbols.

use molframe_core::{BondOrder, Element};
use std::sync::Arc;

/// Broad chemical role of a deposited component.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComponentKind {
    /// Peptide monomer, including modified amino acids.
    AminoAcid,
    /// DNA or RNA monomer.
    Nucleotide,
    /// Carbohydrate monomer.
    Saccharide,
    /// Lipid component.
    Lipid,
    /// Discrete non-polymer component.
    NonPolymer,
    /// Solvent component.
    Solvent,
    /// Monoatomic ion.
    Ion,
    /// A component whose role is not classified by the source.
    #[default]
    Unknown,
}

impl ComponentKind {
    /// Stable integer representation stored in structure annotations.
    #[must_use]
    pub const fn code(self) -> i64 {
        match self {
            Self::Unknown => 0,
            Self::AminoAcid => 1,
            Self::Nucleotide => 2,
            Self::Saccharide => 3,
            Self::Lipid => 4,
            Self::NonPolymer => 5,
            Self::Solvent => 6,
            Self::Ion => 7,
        }
    }

    /// Decodes a structure annotation without interpreting text.
    #[must_use]
    pub const fn from_code(code: i64) -> Option<Self> {
        match code {
            0 => Some(Self::Unknown),
            1 => Some(Self::AminoAcid),
            2 => Some(Self::Nucleotide),
            3 => Some(Self::Saccharide),
            4 => Some(Self::Lipid),
            5 => Some(Self::NonPolymer),
            6 => Some(Self::Solvent),
            7 => Some(Self::Ion),
            _ => None,
        }
    }
}

/// Semantic role of an atom within a polymer component.
///
/// These roles come from a CCD-aware profile or an explicit caller annotation;
/// kernels never reconstruct them from atom-name tables.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PolymerAtomRole(u32);

impl PolymerAtomRole {
    /// No role was assigned.
    pub const UNKNOWN: Self = Self(0);
    /// Protein peptide nitrogen.
    pub const PROTEIN_NITROGEN: Self = Self(1 << 0);
    /// Protein alpha carbon.
    pub const PROTEIN_ALPHA_CARBON: Self = Self(1 << 1);
    /// Protein carbonyl carbon.
    pub const PROTEIN_CARBONYL_CARBON: Self = Self(1 << 2);
    /// Protein carbonyl oxygen.
    pub const PROTEIN_CARBONYL_OXYGEN: Self = Self(1 << 3);
    /// Protein side-chain membership.
    pub const PROTEIN_SIDECHAIN: Self = Self(1 << 4);
    /// Nucleic-acid phosphate.
    pub const NUCLEIC_PHOSPHATE: Self = Self(1 << 5);
    /// Nucleic-acid O5 site.
    pub const NUCLEIC_O5: Self = Self(1 << 6);
    /// Nucleic-acid C5 site.
    pub const NUCLEIC_C5: Self = Self(1 << 7);
    /// Nucleic-acid C4 site.
    pub const NUCLEIC_C4: Self = Self(1 << 8);
    /// Nucleic-acid C3 site.
    pub const NUCLEIC_C3: Self = Self(1 << 9);
    /// Nucleic-acid O3 site.
    pub const NUCLEIC_O3: Self = Self(1 << 10);
    /// Nucleic-acid O4 sugar site.
    pub const NUCLEIC_O4: Self = Self(1 << 11);
    /// Nucleic-acid C1 sugar site.
    pub const NUCLEIC_C1: Self = Self(1 << 12);
    /// Nucleic-acid C2 sugar site.
    pub const NUCLEIC_C2: Self = Self(1 << 13);
    /// Glycosidic nucleobase attachment atom.
    pub const NUCLEIC_GLYCOSIDIC: Self = Self(1 << 14);
    /// Nucleobase atom used to orient the glycosidic torsion/frame.
    pub const NUCLEIC_BASE_REFERENCE: Self = Self(1 << 15);
    /// Other nucleobase membership.
    pub const NUCLEIC_BASE: Self = Self(1 << 16);
    /// Protein beta carbon used to define the first side-chain direction.
    pub const PROTEIN_BETA_CARBON: Self = Self(1 << 17);

    /// Protein backbone membership mask.
    pub const PROTEIN_BACKBONE: Self = Self(
        Self::PROTEIN_NITROGEN.0
            | Self::PROTEIN_ALPHA_CARBON.0
            | Self::PROTEIN_CARBONYL_CARBON.0
            | Self::PROTEIN_CARBONYL_OXYGEN.0,
    );
    /// Nucleic backbone membership mask.
    pub const NUCLEIC_BACKBONE: Self = Self(
        Self::NUCLEIC_PHOSPHATE.0
            | Self::NUCLEIC_O5.0
            | Self::NUCLEIC_C5.0
            | Self::NUCLEIC_C4.0
            | Self::NUCLEIC_C3.0
            | Self::NUCLEIC_O3.0,
    );
    /// Nucleic sugar membership mask.
    pub const NUCLEIC_SUGAR: Self = Self(
        Self::NUCLEIC_C1.0
            | Self::NUCLEIC_C2.0
            | Self::NUCLEIC_C3.0
            | Self::NUCLEIC_C4.0
            | Self::NUCLEIC_O4.0,
    );
    /// Nucleobase membership mask.
    pub const NUCLEIC_BASE_GROUP: Self =
        Self(Self::NUCLEIC_GLYCOSIDIC.0 | Self::NUCLEIC_BASE_REFERENCE.0 | Self::NUCLEIC_BASE.0);
    const ALL: u32 = (1 << 18) - 1;

    /// Stable integer representation stored in atom annotations.
    #[must_use]
    pub const fn code(self) -> i64 {
        self.0 as i64
    }

    /// Decodes a typed atom annotation value.
    #[must_use]
    pub fn from_code(code: i64) -> Option<Self> {
        let Ok(code) = u32::try_from(code) else {
            return None;
        };
        (code & !Self::ALL == 0).then_some(Self(code))
    }

    /// Whether this role set contains any bit in `required`.
    #[must_use]
    pub const fn intersects(self, required: Self) -> bool {
        self.0 & required.0 != 0
    }

    /// Combines independent semantic roles for one atom.
    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

/// Declared tetrahedral configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoConfiguration {
    /// Clockwise priority order.
    R,
    /// Counter-clockwise priority order.
    S,
    /// A mixture or unresolved centre.
    Mixed,
}

/// One expected atom in a chemical component.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComponentAtom {
    /// Component-local atom identifier.
    pub name: Box<str>,
    /// Alternative atom identifier.
    pub alternate_name: Option<Box<str>>,
    /// Chemical element.
    pub element: Element,
    /// Formal charge.
    pub charge: i8,
    /// Aromatic annotation.
    pub aromatic: bool,
    /// Whether the atom leaves during polymerisation.
    pub leaving: bool,
    /// Declared stereochemical configuration.
    pub stereo: Option<StereoConfiguration>,
}

/// One component-local bond.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComponentBond {
    /// First component-local atom identifier.
    pub atom_a: Box<str>,
    /// Second component-local atom identifier.
    pub atom_b: Box<str>,
    /// Chemical order.
    pub order: BondOrder,
    /// Aromatic annotation from the dictionary.
    pub aromatic: bool,
    /// Declared stereochemical configuration.
    pub stereo: Option<StereoConfiguration>,
}

/// A complete deposited-component definition.
#[derive(Clone, Debug, PartialEq)]
pub struct Component {
    /// Dictionary component identifier.
    pub id: Box<str>,
    /// Human-readable name.
    pub name: Box<str>,
    /// Chemical role.
    pub kind: ComponentKind,
    /// Parent standard component for a modification.
    pub parent: Option<Box<str>>,
    /// CCD-declared one-letter polymer code, when it is exactly one ASCII byte.
    pub one_letter_code: Option<u8>,
    /// Formula as deposited by the dictionary.
    pub formula: Option<Box<str>>,
    /// Expected atoms in dictionary order.
    pub atoms: Arc<[ComponentAtom]>,
    /// Internal bonds in dictionary order.
    pub bonds: Arc<[ComponentBond]>,
    /// Ideal coordinates aligned with `atoms`, when complete.
    pub ideal_coordinates: Option<Arc<[[f32; 3]]>>,
    /// Model coordinates aligned with `atoms`, when complete.
    pub model_coordinates: Option<Arc<[[f32; 3]]>>,
}

impl Component {
    /// Finds an expected atom without interpreting its name.
    #[must_use]
    pub fn atom(&self, name: &str) -> Option<&ComponentAtom> {
        self.atoms.iter().find(|atom| atom.name.as_ref() == name)
    }
}
