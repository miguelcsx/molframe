//! Exclusive column storage that survives consecutive batch leases.

use super::StructureBatch;

/// Finished structural columns whose backing allocations can be recycled.
///
/// The wrapper is deliberately not cloneable: a producer keeps exclusive
/// access between leases, so resetting a recycled batch never copies columns.
#[derive(Debug)]
pub struct StructureBatchBuffer {
    batch: StructureBatch,
}

impl StructureBatchBuffer {
    pub(super) fn new(batch: StructureBatch) -> Self {
        Self { batch }
    }

    pub(super) fn batch_mut(&mut self) -> &mut StructureBatch {
        &mut self.batch
    }

    /// Borrows the finished immutable batch.
    #[must_use]
    pub fn batch(&self) -> &StructureBatch {
        &self.batch
    }

    /// Extracts owned columns for a caller that explicitly materialises them.
    #[must_use]
    pub fn into_owned(self) -> StructureBatch {
        self.batch
    }
}
