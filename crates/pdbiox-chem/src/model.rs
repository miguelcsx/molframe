//! Chemical component definitions independent of structure-local symbols.

use pdbiox_core::{BondOrder, Element};
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
