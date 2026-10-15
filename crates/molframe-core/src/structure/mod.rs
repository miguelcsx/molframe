//! The normalised structure, the views over it, and the invariants it holds.
//!
//! A structure is immutable and shared by reference; editing produces a new one.
//! That is what makes handing out a pointer to the coordinates safe, and it is
//! why a view keeps working after the structure it came from has been superseded
//! — it holds its own reference to the snapshot it was made against.

#[cfg(test)]
mod fixture;

mod chain_labels;
mod data;
mod diff;
mod disorder;
mod disorder_checked;
mod edit;
mod ensemble;
mod extensions;
#[cfg(test)]
mod extensions_tests;
mod handle;
mod handle_lookup;
mod label_rewrites;
mod materialize;
mod merge;
mod merge_annotations;
mod overlay;
mod sequence;
mod validate;
mod view;

pub use data::{CoordinateStore, EntryMetadata, Structure, StructureData, UnitCell};
pub use diff::{
    CountDifference, DifferenceError, MetadataDifference, StructureDifference,
    StructureDifferenceOptions, ValueDifference, structure_difference,
};
pub use disorder::{altloc_compatible, compatible_pairs};
pub use edit::CoordinateEditor;
pub use ensemble::EnsembleBlock;
pub use extensions::ExtensionStore;
pub use handle::{AtomRef, ChainRef, ModelRef, ResidueRef};
pub use label_rewrites::LabelRewrites;
pub use overlay::StructureEditor;
pub use sequence::{
    ChainSequenceExt, MissingResidue, ReferenceAlignment, ReferenceSequence,
    SEQUENCE_REFERENCES_EXTENSION, SequenceMapping, SequenceReferences,
};
pub use validate::validate;
pub(crate) use validate::validate_coordinate_edit;
pub use view::StructureView;
