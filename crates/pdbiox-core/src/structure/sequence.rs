//! Observed, canonical and database-reference sequence views.

use super::ChainRef;
use crate::ResidueIndex;
use crate::symbol::SymbolId;

/// Stable extension key for database sequence references and alignments.
pub const SEQUENCE_REFERENCES_EXTENSION: &str = "pdbiox.sequence.references.v1";

/// One database sequence named by `_struct_ref`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceSequence {
    /// Source row identifier.
    pub id: Box<str>,
    /// Entity instantiated by this reference.
    pub entity_id: Box<str>,
    /// Database name, such as `UNP`.
    pub database_name: Option<Box<str>>,
    /// Database code/name.
    pub database_code: Option<Box<str>>,
    /// Stable database accession.
    pub accession: Option<Box<str>>,
    /// Deposited one-letter reference sequence, with CIF whitespace removed.
    pub one_letter_code: Option<Box<str>>,
}

/// One canonical-to-reference interval from `_struct_ref_seq`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceAlignment {
    /// Source alignment identifier.
    pub id: Box<str>,
    /// Referenced `_struct_ref.id`.
    pub reference_id: Box<str>,
    /// Author chain identifiers covered by this alignment.
    pub chain_ids: Box<[Box<str>]>,
    /// Inclusive canonical sequence interval.
    pub canonical: [i32; 2],
    /// Inclusive reference database interval.
    pub reference: [i32; 2],
}

/// Typed reference sequence metadata attached to a structure.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SequenceReferences {
    /// Database sequences in source order.
    pub sequences: Vec<ReferenceSequence>,
    /// Canonical-to-reference alignments in source order.
    pub alignments: Vec<ReferenceAlignment>,
}

/// One observed residue's correspondence across sequence namespaces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SequenceMapping<'a> {
    /// Observed structure residue.
    pub residue: ResidueIndex,
    /// Canonical entity position from `label_seq_id`.
    pub canonical_position: Option<i32>,
    /// Reference database position when an alignment covers it.
    pub reference_position: Option<i32>,
    /// Reference identifier supplying that position.
    pub reference_id: Option<&'a str>,
}

/// One canonical residue absent from the observed structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MissingResidue {
    /// One-based canonical sequence position.
    pub canonical_position: u32,
    /// Expected component identifier.
    pub component: SymbolId,
}

/// Sequence views carried by a hierarchy chain handle.
pub trait ChainSequenceExt<'a> {
    /// Component identifiers actually modelled, in topology order.
    fn observed_sequence(self) -> Vec<SymbolId>;
    /// Entity sequence that should have been present.
    fn canonical_sequence(self) -> &'a [SymbolId];
    /// Database sequences associated with the chain's entity.
    fn reference_sequences(self) -> Vec<&'a ReferenceSequence>;
    /// Observed-to-canonical-to-reference position mapping.
    fn sequence_mapping(self) -> Vec<SequenceMapping<'a>>;
    /// Canonical positions not represented by an observed residue.
    fn missing_residues(self) -> Vec<MissingResidue>;
}

impl<'a> ChainSequenceExt<'a> for ChainRef<'a> {
    fn observed_sequence(self) -> Vec<SymbolId> {
        self.residues()
            .filter_map(super::ResidueRef::label_comp_id)
            .collect()
    }

    fn canonical_sequence(self) -> &'a [SymbolId] {
        let Some(entity) = self.entity() else {
            return &[];
        };
        self.data.topology.entities.canonical_sequence(entity)
    }

    fn reference_sequences(self) -> Vec<&'a ReferenceSequence> {
        let Some(entity) = self.entity() else {
            return Vec::new();
        };
        let Some(entity_id) = self
            .data
            .topology
            .entities
            .id(entity)
            .and_then(|id| self.data.dictionary.resolve(id))
        else {
            return Vec::new();
        };
        match self.references() {
            Some(references) => references
                .sequences
                .iter()
                .filter(|sequence| sequence.entity_id.as_ref() == entity_id)
                .collect(),
            None => Vec::new(),
        }
    }

    fn sequence_mapping(self) -> Vec<SequenceMapping<'a>> {
        let alignment = self.alignment();
        self.residues()
            .map(|residue| {
                let canonical_position = residue.label_seq_id();
                let (reference_position, reference_id) = alignment
                    .and_then(|alignment| {
                        let position = canonical_position?;
                        if position < alignment.canonical[0] || position > alignment.canonical[1] {
                            return None;
                        }
                        Some((
                            alignment.reference[0] + position - alignment.canonical[0],
                            alignment.reference_id.as_ref(),
                        ))
                    })
                    .map_or((None, None), |(position, id)| (Some(position), Some(id)));
                SequenceMapping {
                    residue: residue.index(),
                    canonical_position,
                    reference_position,
                    reference_id,
                }
            })
            .collect()
    }

    fn missing_residues(self) -> Vec<MissingResidue> {
        let observed = self
            .residues()
            .filter_map(super::ResidueRef::label_seq_id)
            .collect::<std::collections::BTreeSet<_>>();
        self.canonical_sequence()
            .iter()
            .enumerate()
            .filter_map(|(position, component)| {
                let canonical_position = u32::try_from(position).ok()?.checked_add(1)?;
                let observed_position = i32::try_from(canonical_position).ok()?;
                (!observed.contains(&observed_position)).then_some(MissingResidue {
                    canonical_position,
                    component: *component,
                })
            })
            .collect()
    }
}

impl<'a> ChainRef<'a> {
    fn references(self) -> Option<&'a SequenceReferences> {
        self.data.extensions.get(SEQUENCE_REFERENCES_EXTENSION)
    }

    fn alignment(self) -> Option<&'a ReferenceAlignment> {
        let label = self.label();
        let auth = self.auth_label();
        self.references()?.alignments.iter().find(|alignment| {
            alignment
                .chain_ids
                .iter()
                .any(|chain| label == Some(chain.as_ref()) || auth == Some(chain.as_ref()))
        })
    }
}
