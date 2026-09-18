//! Explicit expansion of a lazy biological assembly.

use crate::AssemblyView;
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation, AtomAnnotations};
use molframe_core::bond::{BondRecord, BondTableBuilder};
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::hashing::{IdentityHashMap, IdentityHashSet};
use molframe_core::optional::OptionalSymbol;
use molframe_core::structure::CoordinateStore;
use molframe_core::symbol::SymbolId;
use molframe_core::topology::{ChainRecord, ResidueRecord};
use molframe_core::{
    AtomIndex, ChainIndex, Code, CoordinateBlock, CoordinateGeneration, Diagnostic, ModelIndex,
    ResidueIndex, Structure, StructureData,
};
use molframe_geom::Rigid;
use std::collections::BTreeMap;
use std::ops::Range;

/// Per-atom annotation holding the generated chain's stable instance identifier.
pub const INSTANCE_ID_ANNOTATION: &str = "molframe.instance_id";

/// Describes one materialized copy of a source-chain atom span.
#[derive(Clone, Debug)]
struct CopySpan {
    source_chain: ChainIndex,
    transform: Rigid,
    source_atoms: Range<u32>,
    output_start: u32,
}

/// Groups copy indices by source chain for one identical rigid transform.
type CopyGroup = BTreeMap<u32, Vec<usize>>;

/// Exact bit representation of a rigid transform used for deterministic grouping.
type TransformKey = [u64; 12];

include!("materialize/impl.rs");
include!("materialize/helpers.rs");

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;
