//! MMTF-only metadata that the normalized structure model cannot express.

/// Stable extension key for metadata required by lossless MMTF writing.
pub const MMTF_METADATA_EXTENSION: &str = "molframe.mmtf.metadata.v1";

/// Chemistry fields for one MMTF group dictionary entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MmtfGroupMetadata {
    /// Component identifier.
    pub name: Box<str>,
    /// Atom names used to verify that this dictionary entry still applies.
    pub atom_names: Vec<Box<str>>,
    /// Deposited element symbols, or absence of the complete optional list.
    pub elements: Option<Vec<Box<str>>>,
    /// Deposited one-letter code.
    pub single_letter_code: Box<str>,
    /// Deposited `PDBx` chemical component type.
    pub chem_comp_type: Box<str>,
}

/// Metadata for one MMTF entity entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MmtfEntityMetadata {
    /// Description, including an explicitly empty description.
    pub description: Box<str>,
    /// Deposited MMTF entity type.
    pub kind: Box<str>,
    /// Deposited one-letter entity sequence.
    pub sequence: Box<str>,
}

/// Source fields retained for faithful MMTF re-encoding.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MmtfMetadata {
    /// Deposited space-group name.
    pub space_group: Option<Box<str>>,
    /// Group dictionary chemistry in source order.
    pub groups: Vec<MmtfGroupMetadata>,
    /// Entity metadata in source order.
    pub entities: Vec<MmtfEntityMetadata>,
    /// Optional columns present in the source document.
    pub optional_fields: std::collections::BTreeSet<MmtfOptionalField>,
}

/// One optional top-level MMTF column.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum MmtfOptionalField {
    /// Per-atom B factors.
    BFactor,
    /// Per-atom occupancies.
    Occupancy,
    /// Per-atom source identifiers.
    AtomId,
    /// Per-atom alternate locations.
    AltLoc,
    /// Per-group insertion codes.
    InsCode,
    /// Per-group canonical sequence indices.
    SequenceIndex,
    /// Per-chain author names.
    ChainName,
    /// Entity metadata and chain membership.
    EntityList,
}
