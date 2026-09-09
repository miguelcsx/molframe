//! Explicit expansion of a lazy biological assembly.

use crate::AssemblyView;
use pdbiox_core::annotation::{AnnotationColumn, AtomAnnotation, AtomAnnotations};
use pdbiox_core::bond::{BondRecord, BondTableBuilder};
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::hashing::{IdentityHashMap, IdentityHashSet};
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::structure::CoordinateStore;
use pdbiox_core::symbol::SymbolId;
use pdbiox_core::topology::{ChainRecord, ResidueRecord};
use pdbiox_core::{
    AtomIndex, ChainIndex, Code, CoordinateBlock, CoordinateGeneration, Diagnostic, ModelIndex,
    ResidueIndex, Structure, StructureData,
};
use pdbiox_geom::Rigid;
use std::collections::BTreeMap;
use std::ops::Range;

/// Per-atom annotation holding the generated chain's stable instance identifier.
pub const INSTANCE_ID_ANNOTATION: &str = "pdbiox.instance_id";

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
