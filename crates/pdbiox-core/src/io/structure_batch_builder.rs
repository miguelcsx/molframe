//! Allocation-aware construction and reuse of structure batches.

use super::{BatchContinuity, StructureAtomRecord, StructureBatch, StructureBatchBuffer};
use crate::Diagnostic;
use crate::column::Presence;
use crate::optional::OptionalSymbol;
use crate::provider::{ChunkDescriptor, ChunkId, DatasetId, LogicalRow, ProviderError};
use crate::symbol::{AltId, DictionaryFull, Interner};

/// Incremental builder that interns each distinct identifier once per batch.
#[derive(Debug)]
pub struct StructureBatchBuilder {
    buffer: StructureBatchBuffer,
}

impl StructureBatchBuilder {
    /// Creates an empty batch with a bounded row reservation.
    ///
    /// # Errors
    ///
    /// Rejects an unrepresentable descriptor range.
    pub fn new(
        dataset: DatasetId,
        chunk: ChunkId,
        logical_start: LogicalRow,
        capacity: u32,
    ) -> Result<Self, ProviderError> {
        let descriptor = ChunkDescriptor::new(dataset, chunk, logical_start, 0)?;
        let capacity = capacity as usize;
        Ok(Self {
            buffer: StructureBatchBuffer::new(StructureBatch {
                descriptor,
                continuity: BatchContinuity::default(),
                dictionary: Interner::new(),
                model: Vec::with_capacity(capacity),
                chain: Vec::with_capacity(capacity),
                component: Vec::with_capacity(capacity),
                sequence: Vec::with_capacity(capacity),
                insertion: Vec::with_capacity(capacity),
                atom: Vec::with_capacity(capacity),
                alternate: Vec::with_capacity(capacity),
                element: Vec::with_capacity(capacity),
                position: Vec::with_capacity(capacity),
                coordinate_presence: Vec::with_capacity(capacity),
                occupancy: Vec::with_capacity(capacity),
                occupancy_presence: Vec::with_capacity(capacity),
                b_factor: Vec::with_capacity(capacity),
                b_factor_presence: Vec::with_capacity(capacity),
                formal_charge: Vec::with_capacity(capacity),
                formal_charge_presence: Vec::with_capacity(capacity),
                atom_site_id: Vec::with_capacity(capacity),
                heterogen: Vec::with_capacity(capacity),
                diagnostics: Vec::new(),
            }),
        })
    }

    /// Resets an exclusively owned batch while preserving column capacity.
    ///
    /// The caller must first use [`StructureBatch::can_reuse`] to ensure the
    /// requested row extent needs no allocation.
    ///
    /// # Errors
    ///
    /// Rejects an unrepresentable descriptor range.
    pub(super) fn reuse(
        mut buffer: StructureBatchBuffer,
        dataset: DatasetId,
        chunk: ChunkId,
        logical_start: LogicalRow,
    ) -> Result<Self, ProviderError> {
        let batch = buffer.batch_mut();
        batch.descriptor = ChunkDescriptor::new(dataset, chunk, logical_start, 0)?;
        batch.continuity = BatchContinuity::default();
        batch.dictionary.clear();
        batch.model.clear();
        batch.chain.clear();
        batch.component.clear();
        batch.sequence.clear();
        batch.insertion.clear();
        batch.atom.clear();
        batch.alternate.clear();
        batch.element.clear();
        batch.position.clear();
        batch.coordinate_presence.clear();
        batch.occupancy.clear();
        batch.occupancy_presence.clear();
        batch.b_factor.clear();
        batch.b_factor_presence.clear();
        batch.formal_charge.clear();
        batch.formal_charge_presence.clear();
        batch.atom_site_id.clear();
        batch.heterogen.clear();
        batch.diagnostics.clear();
        Ok(Self { buffer })
    }

    /// Sets hierarchy continuity discovered by the source.
    pub fn set_continuity(&mut self, continuity: BatchContinuity) {
        self.buffer.batch_mut().continuity = continuity;
    }

    /// Appends one row without allocating per atom.
    ///
    /// # Errors
    ///
    /// Returns [`DictionaryFull`] if the batch-local dictionary is exhausted.
    pub fn push(&mut self, record: StructureAtomRecord<'_>) -> Result<(), DictionaryFull> {
        let batch = self.buffer.batch_mut();
        let chain = batch.dictionary.intern(record.chain)?;
        let component = batch.dictionary.intern(record.component)?;
        let atom = batch.dictionary.intern(record.atom)?;
        let insertion = optional_symbol(&mut batch.dictionary, record.insertion)?;
        let alternate = optional_alt(&mut batch.dictionary, record.alternate)?;
        batch.model.push(record.model);
        batch.chain.push(chain);
        batch.component.push(component);
        batch.sequence.push(record.sequence);
        batch.insertion.push(insertion);
        batch.atom.push(atom);
        batch.alternate.push(alternate);
        batch.element.push(record.element.atomic_number());
        if let Some(position) = record.position {
            batch.position.push(position);
            batch.coordinate_presence.push(Presence::Present);
        } else {
            batch.position.push([f32::NAN; 3]);
            batch.coordinate_presence.push(Presence::Unknown);
        }
        batch.occupancy.push(record.occupancy.0);
        batch.occupancy_presence.push(record.occupancy.1);
        batch.b_factor.push(record.b_factor.0);
        batch.b_factor_presence.push(record.b_factor.1);
        batch.formal_charge.push(record.formal_charge.0);
        batch.formal_charge_presence.push(record.formal_charge.1);
        batch.atom_site_id.push(record.atom_site_id);
        batch.heterogen.push(u8::from(record.heterogen));
        Ok(())
    }

    /// Appends a logical-row diagnostic.
    pub fn push_diagnostic(&mut self, diagnostic: Diagnostic) {
        self.buffer.batch_mut().diagnostics.push(diagnostic);
    }

    /// Finalises the descriptor using the actual local extent.
    ///
    /// # Errors
    ///
    /// Rejects more than `u32::MAX` rows or logical identity overflow.
    pub fn finish(mut self) -> Result<StructureBatchBuffer, ProviderError> {
        let batch = self.buffer.batch_mut();
        let rows =
            u32::try_from(batch.model.len()).map_err(|_| ProviderError::IdentityOverflow {
                dataset: batch.descriptor.dataset(),
                field: "batch rows",
            })?;
        batch.descriptor = ChunkDescriptor::new(
            batch.descriptor.dataset(),
            batch.descriptor.chunk(),
            batch.descriptor.logical_start(),
            rows,
        )?;
        Ok(self.buffer)
    }
}

fn optional_symbol(
    dictionary: &mut Interner,
    text: &str,
) -> Result<OptionalSymbol, DictionaryFull> {
    if text.is_empty() {
        return Ok(OptionalSymbol::NONE);
    }
    dictionary.intern(text).map(OptionalSymbol::some)
}

fn optional_alt(dictionary: &mut Interner, text: &str) -> Result<AltId, DictionaryFull> {
    if text.is_empty() {
        return Ok(AltId::BLANK);
    }
    let symbol = dictionary.intern(text)?;
    AltId::labelled(symbol).ok_or(DictionaryFull)
}
