//! Complete Tripos metadata retained alongside the shared molecular graph.

use crate::Molecule;

/// Metadata aligned with one MOL2 atom.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2AtomMetadata {
    /// Positive source atom identifier.
    pub id: usize,
    /// Tripos atom name.
    pub name: Box<str>,
    /// SYBYL atom type.
    pub atom_type: Box<str>,
    /// Optional substructure identifier.
    pub substructure_id: Option<usize>,
    /// Optional substructure name.
    pub substructure_name: Option<Box<str>>,
    /// Optional partial charge.
    pub charge: Option<f64>,
    /// Optional status bit expression.
    pub status_bits: Option<Box<str>>,
}

/// Metadata aligned with one MOL2 bond.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mol2BondMetadata {
    /// Positive source bond identifier.
    pub id: usize,
    /// Tripos bond type, including `ar`, `am`, `du`, `un`, or `nc`.
    pub bond_type: Box<str>,
    /// Optional status bit expression.
    pub status_bits: Option<Box<str>>,
}

/// An additional Tripos section retained verbatim by line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mol2Section {
    /// Section name without the `@<TRIPOS>` prefix.
    pub name: Box<str>,
    /// Section body in source order.
    pub lines: Vec<Box<str>>,
}

/// One complete MOL2 molecule record.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2Record {
    /// Molecule name.
    pub name: Box<str>,
    /// Tripos molecule type.
    pub molecule_type: Box<str>,
    /// Tripos charge type.
    pub charge_type: Box<str>,
    /// Counts following atom and bond counts, preserved at their source arity.
    pub additional_counts: Vec<usize>,
    /// Optional molecule status bits.
    pub status_bits: Option<Box<str>>,
    /// Optional molecule comment.
    pub comment: Option<Box<str>>,
    /// Shared molecular graph.
    pub molecule: Molecule,
    /// Atom metadata aligned with the graph.
    pub atom_metadata: Vec<Mol2AtomMetadata>,
    /// Bond metadata aligned with the graph.
    pub bond_metadata: Vec<Mol2BondMetadata>,
    /// All other sections in source order.
    pub extra_sections: Vec<Mol2Section>,
}

/// Malformed or unrepresentable MOL2 data.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Mol2Error {
    /// Required syntax or a number is malformed.
    Malformed,
    /// Declared counts disagree with parsed sections.
    CountMismatch,
    /// Parallel metadata does not match the graph.
    MetadataLength,
    /// An identifier is zero, duplicated, or unresolved.
    InvalidIdentifier,
    /// A field cannot be serialized without corrupting MOL2 syntax.
    Unrepresentable,
}

impl std::fmt::Display for Mol2Error {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let message = match self {
            Self::Malformed => "malformed MOL2 record",
            Self::CountMismatch => "MOL2 declared counts do not match its sections",
            Self::MetadataLength => "MOL2 metadata length does not match atoms or bonds",
            Self::InvalidIdentifier => "MOL2 identifier is zero, duplicated, or unresolved",
            Self::Unrepresentable => "value cannot be represented by MOL2",
        };
        formatter.write_str(message)
    }
}

impl std::error::Error for Mol2Error {}
