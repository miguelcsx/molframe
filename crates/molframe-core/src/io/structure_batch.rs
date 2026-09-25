//! Bounded columnar rows shared by structural format readers.

use crate::column::Presence;
use crate::execution::Batch;
use crate::optional::OptionalSymbol;
use crate::provider::{ChunkDescriptor, ProviderError};
use crate::symbol::{AltId, DictionaryFull, Interner, SymbolId};
use crate::{Diagnostic, Element};
use std::fmt;

/// Failure shared by bounded structural batch readers and collectors.
#[derive(Debug)]
pub enum StructureBatchError {
    /// A source, syntax or semantic diagnostic.
    Diagnostic(Diagnostic),
    /// A single record cannot fit the caller's row or byte demand.
    DemandTooSmall {
        /// Minimum retained bytes required.
        required: usize,
        /// Retained bytes offered by the demand.
        available: usize,
    },
    /// A single indivisible record cannot fit the execution budget.
    RecordExceedsBudget {
        /// Minimum bytes required for the indivisible record.
        required: usize,
        /// Bytes available to the source stage.
        available: usize,
    },
    /// Stable dataset, chunk or logical-row arithmetic overflowed.
    Identity(ProviderError),
    /// Execution memory was exhausted after source advancement was prevented.
    Memory(crate::MemoryBudgetError),
    /// A batch-local identifier dictionary was exhausted.
    DictionaryFull,
}

impl fmt::Display for StructureBatchError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Diagnostic(error) => error.fmt(formatter),
            Self::DemandTooSmall {
                required,
                available,
            } => write!(
                formatter,
                "record requires {required} bytes but demand permits {available}"
            ),
            Self::RecordExceedsBudget {
                required,
                available,
            } => write!(
                formatter,
                "record requires {required} bytes but budget permits {available}"
            ),
            Self::Identity(error) => error.fmt(formatter),
            Self::Memory(error) => error.fmt(formatter),
            Self::DictionaryFull => formatter.write_str("batch identifier dictionary is full"),
        }
    }
}

impl std::error::Error for StructureBatchError {}

impl From<Diagnostic> for StructureBatchError {
    fn from(error: Diagnostic) -> Self {
        Self::Diagnostic(error)
    }
}

impl From<ProviderError> for StructureBatchError {
    fn from(error: ProviderError) -> Self {
        Self::Identity(error)
    }
}

impl From<DictionaryFull> for StructureBatchError {
    fn from(_: DictionaryFull) -> Self {
        Self::DictionaryFull
    }
}

impl From<crate::MemoryBudgetError> for StructureBatchError {
    fn from(error: crate::MemoryBudgetError) -> Self {
        Self::Memory(error)
    }
}

/// Deepest hierarchy level that continues across one batch boundary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ContinuityLevel {
    /// No hierarchy row crosses the boundary.
    #[default]
    None,
    /// The model continues but its chain does not.
    Model,
    /// The model and chain continue but the residue does not.
    Chain,
    /// The model, chain and residue all continue.
    Residue,
}

/// Whether hierarchy rows continue across adjacent batches.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BatchContinuity {
    /// Deepest hierarchy row begun in the preceding batch.
    pub before: ContinuityLevel,
    /// Deepest hierarchy row completed in the following batch.
    pub after: ContinuityLevel,
}

/// One atom row accepted by [`StructureBatchBuilder`](crate::io::StructureBatchBuilder).
#[derive(Clone, Copy, Debug)]
pub struct StructureAtomRecord<'a> {
    /// Deposited model number.
    pub model: i32,
    /// Normalised chain identifier.
    pub chain: &'a str,
    /// Component identifier.
    pub component: &'a str,
    /// Deposited residue sequence number, or `i32::MIN` when absent.
    pub sequence: i32,
    /// Insertion code, or an empty string when absent.
    pub insertion: &'a str,
    /// Atom identifier.
    pub atom: &'a str,
    /// Alternate-location identifier, or an empty string when absent.
    pub alternate: &'a str,
    /// Chemical element.
    pub element: Element,
    /// Cartesian position, or `None` when unrecorded.
    pub position: Option<[f32; 3]>,
    /// Occupancy and validity.
    pub occupancy: (f32, Presence),
    /// B factor and validity.
    pub b_factor: (f32, Presence),
    /// Formal charge and validity.
    pub formal_charge: (i8, Presence),
    /// Deposited atom-site identifier, with zero reserved for absence.
    pub atom_site_id: u32,
    /// Whether the source marked the residue as a heterogen.
    pub heterogen: bool,
}

/// Bounded structure rows in a format-neutral layout.
#[derive(Clone, Debug)]
pub struct StructureBatch {
    pub(super) descriptor: ChunkDescriptor,
    pub(super) continuity: BatchContinuity,
    pub(super) dictionary: Interner,
    pub(super) model: Vec<i32>,
    pub(super) chain: Vec<SymbolId>,
    pub(super) component: Vec<SymbolId>,
    pub(super) sequence: Vec<i32>,
    pub(super) insertion: Vec<OptionalSymbol>,
    pub(super) atom: Vec<SymbolId>,
    pub(super) alternate: Vec<AltId>,
    pub(super) element: Vec<u8>,
    pub(super) position: Vec<[f32; 3]>,
    pub(super) coordinate_presence: Vec<Presence>,
    pub(super) occupancy: Vec<f32>,
    pub(super) occupancy_presence: Vec<Presence>,
    pub(super) b_factor: Vec<f32>,
    pub(super) b_factor_presence: Vec<Presence>,
    pub(super) formal_charge: Vec<i8>,
    pub(super) formal_charge_presence: Vec<Presence>,
    pub(super) atom_site_id: Vec<u32>,
    pub(super) heterogen: Vec<u8>,
    pub(super) diagnostics: Vec<Diagnostic>,
}

impl StructureBatch {
    /// Reports whether this allocation can serve another batch without
    /// growing any hot row column or exceeding the next byte demand.
    #[must_use]
    pub fn can_reuse(&self, rows: u32, max_retained_bytes: usize) -> bool {
        let rows = rows as usize;
        self.model.capacity() >= rows && self.retained_bytes() <= max_retained_bytes
    }

    /// Stable global identity and local row extent.
    #[must_use]
    pub const fn descriptor(&self) -> ChunkDescriptor {
        self.descriptor
    }

    /// Hierarchy continuation flags.
    #[must_use]
    pub const fn continuity(&self) -> BatchContinuity {
        self.continuity
    }

    /// Batch-local identifier dictionary.
    #[must_use]
    pub const fn dictionary(&self) -> &Interner {
        &self.dictionary
    }

    /// Deposited model numbers.
    #[must_use]
    pub fn models(&self) -> &[i32] {
        &self.model
    }

    /// Chain identifiers.
    #[must_use]
    pub fn chains(&self) -> &[SymbolId] {
        &self.chain
    }

    /// Component identifiers.
    #[must_use]
    pub fn components(&self) -> &[SymbolId] {
        &self.component
    }

    /// Deposited sequence numbers, with `i32::MIN` reserved for absence.
    #[must_use]
    pub fn sequences(&self) -> &[i32] {
        &self.sequence
    }

    /// Insertion codes.
    #[must_use]
    pub fn insertions(&self) -> &[OptionalSymbol] {
        &self.insertion
    }

    /// Atom identifiers.
    #[must_use]
    pub fn atoms(&self) -> &[SymbolId] {
        &self.atom
    }

    /// Alternate-location identifiers.
    #[must_use]
    pub fn alternates(&self) -> &[AltId] {
        &self.alternate
    }

    /// Atomic numbers.
    #[must_use]
    pub fn elements(&self) -> &[u8] {
        &self.element
    }

    /// Cartesian positions; missing rows contain `NaN`.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        &self.position
    }

    /// Position validity.
    #[must_use]
    pub fn coordinate_presence(&self) -> &[Presence] {
        &self.coordinate_presence
    }

    /// Occupancy values and validity.
    #[must_use]
    pub fn occupancies(&self) -> (&[f32], &[Presence]) {
        (&self.occupancy, &self.occupancy_presence)
    }

    /// B factors and validity.
    #[must_use]
    pub fn b_factors(&self) -> (&[f32], &[Presence]) {
        (&self.b_factor, &self.b_factor_presence)
    }

    /// Formal charges and validity.
    #[must_use]
    pub fn formal_charges(&self) -> (&[i8], &[Presence]) {
        (&self.formal_charge, &self.formal_charge_presence)
    }

    /// Deposited atom-site identifiers.
    #[must_use]
    pub fn atom_site_ids(&self) -> &[u32] {
        &self.atom_site_id
    }

    /// Per-row heterogen sentinels.
    #[must_use]
    pub fn heterogens(&self) -> &[u8] {
        &self.heterogen
    }

    /// Diagnostics emitted while producing these logical rows.
    #[must_use]
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }
}

impl Batch for StructureBatch {
    fn rows(&self) -> usize {
        self.model.len()
    }

    fn retained_bytes(&self) -> usize {
        let mut bytes = self.dictionary.retained_bytes();
        bytes = bytes.saturating_add(self.model.capacity_bytes());
        bytes = bytes.saturating_add(self.chain.capacity_bytes());
        bytes = bytes.saturating_add(self.component.capacity_bytes());
        bytes = bytes.saturating_add(self.sequence.capacity_bytes());
        bytes = bytes.saturating_add(self.insertion.capacity_bytes());
        bytes = bytes.saturating_add(self.atom.capacity_bytes());
        bytes = bytes.saturating_add(self.alternate.capacity_bytes());
        bytes = bytes.saturating_add(self.element.capacity_bytes());
        bytes = bytes.saturating_add(self.position.capacity_bytes());
        bytes = bytes.saturating_add(self.coordinate_presence.capacity_bytes());
        bytes = bytes.saturating_add(self.occupancy.capacity_bytes());
        bytes = bytes.saturating_add(self.occupancy_presence.capacity_bytes());
        bytes = bytes.saturating_add(self.b_factor.capacity_bytes());
        bytes = bytes.saturating_add(self.b_factor_presence.capacity_bytes());
        bytes = bytes.saturating_add(self.formal_charge.capacity_bytes());
        bytes = bytes.saturating_add(self.formal_charge_presence.capacity_bytes());
        bytes = bytes.saturating_add(self.atom_site_id.capacity_bytes());
        bytes = bytes.saturating_add(self.heterogen.capacity_bytes());
        bytes = bytes.saturating_add(self.diagnostics.capacity_bytes());
        bytes.saturating_add(
            self.diagnostics
                .iter()
                .map(Diagnostic::retained_bytes)
                .sum::<usize>(),
        )
    }
}

trait VecCapacityBytes {
    fn capacity_bytes(&self) -> usize;
}

impl<T> VecCapacityBytes for Vec<T> {
    fn capacity_bytes(&self) -> usize {
        self.capacity().saturating_mul(std::mem::size_of::<T>())
    }
}
